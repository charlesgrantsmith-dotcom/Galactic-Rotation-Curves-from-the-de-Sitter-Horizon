#!/usr/bin/env python3
"""
Robustness check: rerun SPARC fits with different error floors.
"""
import os, sys, warnings
import numpy as np
from scipy.optimize import minimize
from collections import OrderedDict

warnings.filterwarnings('ignore')

C_MS = 2.998e8
KPC_M = 3.086e19
G_SI = 6.674e-11
MSUN = 1.989e30

def au_from_H0(H0):
    H0_si = H0 * 1e3 / 3.086e22
    return C_MS * H0_si / (2 * np.pi)

def g_holographic(gN, au):
    gN = np.asarray(gN, dtype=float)
    x = np.sqrt(np.abs(gN) / au)
    with np.errstate(over='ignore', invalid='ignore'):
        denom = 1.0 - np.exp(-x)
    safe = x > 1e-8
    return np.where(safe, gN / np.where(safe, denom, 1.0),
                    np.sign(gN + 1e-30) * np.sqrt(np.abs(gN) * au))

def g_mond_simple(gN, a0):
    gN = np.asarray(gN, dtype=float)
    return 0.5 * (gN + np.sqrt(gN**2 + 4*np.abs(gN)*a0))

def g_mond_standard(gN, a0):
    gN = np.asarray(gN, dtype=float)
    gN_abs = np.abs(gN)
    g = np.sqrt(gN_abs**2 + gN_abs * a0)
    for _ in range(30):
        x = g / a0
        mu = x / np.sqrt(1 + x**2)
        f = g * mu - gN_abs
        mu_p = 1.0 / (1 + x**2)**1.5
        df = mu + g * mu_p / a0
        g = np.maximum(g - f / (df + 1e-30), 1e-30)
    return g * np.sign(gN + 1e-30)

def load_mass_models(filepath):
    galaxies = OrderedDict()
    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if (line.startswith('Title') or line.startswith('Author') or
                line.startswith('Table') or line.startswith('=') or
                line.startswith('Byte') or line.startswith('-') or
                line.startswith(' ') and ('Bytes' in line or 'Format' in line) or
                line.startswith('Note') or line.strip() == '' or
                'Label' in line or 'Explanations' in line or
                'Units' in line or 'uncertainty' in line or
                'factor' in line or 'solMass' in line or
                'and Accurate' in line):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            try:
                name = parts[0]
                R = float(parts[2]); Vobs = float(parts[3]); eV = float(parts[4])
                Vgas = float(parts[5]); Vdisk = float(parts[6])
                Vbul = float(parts[7]) if len(parts) > 7 else 0.0
            except (ValueError, IndexError):
                continue
            if name not in galaxies:
                galaxies[name] = {'R':[], 'Vobs':[], 'eV':[], 'Vgas':[], 'Vdisk':[], 'Vbul':[]}
            galaxies[name]['R'].append(R)
            galaxies[name]['Vobs'].append(Vobs)
            galaxies[name]['eV'].append(eV)
            galaxies[name]['Vgas'].append(Vgas)
            galaxies[name]['Vdisk'].append(Vdisk)
            galaxies[name]['Vbul'].append(Vbul)
    for name in galaxies:
        for key in galaxies[name]:
            galaxies[name][key] = np.array(galaxies[name][key])
    return galaxies

def compute_gN(R, Vgas, Vdisk, Vbul, Yd, Yb):
    V2 = Yd * Vdisk * np.abs(Vdisk) + Yb * Vbul * np.abs(Vbul) + Vgas * np.abs(Vgas)
    return V2 * 1e6 / (R * KPC_M)

def fit_galaxy(gal, model_func, accel, error_floor_pct, error_floor_abs):
    R, Vobs, eV = gal['R'], gal['Vobs'], gal['eV']
    Vgas, Vdisk, Vbul = gal['Vgas'], gal['Vdisk'], gal['Vbul']

    # Apply error floor
    if error_floor_pct is None and error_floor_abs is None:
        eV_use = np.maximum(eV, 0.01)  # just prevent division by zero
    else:
        floor_abs = error_floor_abs if error_floor_abs else 0
        floor_pct = error_floor_pct * np.abs(Vobs) / 100.0 if error_floor_pct else 0
        eV_use = np.maximum(eV, np.maximum(floor_abs, floor_pct))

    has_bulge = np.any(np.abs(Vbul) > 0.5)

    def chi2(params):
        Yd = params[0]
        Yb = params[1] if has_bulge else 0.0
        if Yd < 0.01 or Yd > 10: return 1e10
        if has_bulge and (Yb < 0.01 or Yb > 10): return 1e10
        gN = compute_gN(R, Vgas, Vdisk, Vbul, Yd, Yb)
        g_eff = model_func(gN, accel)
        Vpred = np.sqrt(np.abs(g_eff) * R * KPC_M) / 1e3
        return np.sum(((Vobs - Vpred) / eV_use)**2)

    starts = [[0.3], [0.5], [0.8], [1.2]]
    bounds = [(0.01, 10.0)]
    if has_bulge:
        starts = [[0.3,0.3],[0.5,0.5],[0.5,0.8],[0.8,0.5],[1.0,1.0]]
        bounds = [(0.01,10.0),(0.01,10.0)]

    best = None
    for x0 in starts:
        try:
            res = minimize(chi2, x0, method='L-BFGS-B', bounds=bounds)
            if best is None or res.fun < best.fun: best = res
        except: pass

    if best is None: return None
    dof = max(len(Vobs) - (2 if has_bulge else 1), 1)
    return {'chi2': best.fun, 'dof': dof, 'chi2r': best.fun/dof, 'Yd': best.x[0],
            'Yb': best.x[1] if has_bulge else 0.0}

def run_scenario(galaxies, good_names, models, floor_label, floor_pct, floor_abs):
    results = {m: [] for m in models}
    for gname in good_names:
        gal = galaxies[gname]
        for mname, (mfunc, accel) in models.items():
            res = fit_galaxy(gal, mfunc, accel, floor_pct, floor_abs)
            if res:
                res['name'] = gname
                results[mname].append(res)
    return results

def print_results(label, results, models):
    print(f"\n{'='*75}")
    print(f"  {label}")
    print(f"{'='*75}")
    print(f"{'Model':<18} {'Med χ²/ν':>9} {'Mean':>7} {'<1':>6} {'<2':>6} {'<5':>6} "
          f"{'Med Υd':>7} {'N':>4}")
    print('-'*70)

    stats = {}
    for mname in models:
        rl = results[mname]
        if not rl: continue
        c = np.array([r['chi2r'] for r in rl])
        y = np.array([r['Yd'] for r in rl])
        stats[mname] = {
            'med': np.median(c), 'mean': np.mean(c),
            'f1': np.mean(c<1)*100, 'f2': np.mean(c<2)*100, 'f5': np.mean(c<5)*100,
            'med_Yd': np.median(y), 'n': len(rl)
        }
        s = stats[mname]
        print(f"{mname:<18} {s['med']:>9.2f} {s['mean']:>7.2f} "
              f"{s['f1']:>5.1f}% {s['f2']:>5.1f}% {s['f5']:>5.1f}% "
              f"{s['med_Yd']:>7.3f} {s['n']:>4}")

    # Head to head
    h = {r['name']: r['chi2r'] for r in results.get('Holographic', [])}
    m = {r['name']: r['chi2r'] for r in results.get('MOND_simple', [])}
    common = sorted(set(h.keys()) & set(m.keys()))
    if common:
        delta = np.array([h[n] - m[n] for n in common])
        wins = np.sum(delta < 0)
        losses = np.sum(delta > 0)
        print(f"\n  Head-to-head vs MOND simple: Holographic wins {wins}/{len(common)} "
              f"({wins/len(common)*100:.1f}%), median Δχ²/ν = {np.median(delta):.3f}")

    return stats

def main():
    H0 = 67.4
    a_u = au_from_H0(H0)
    a0 = 1.2e-10

    print(f"a_u = {a_u:.4e} m/s²  |  a0 = {a0:.4e} m/s²  |  ratio = {a_u/a0:.4f}")

    galaxies = load_mass_models('/home/claude/SPARC/MassModels_Lelli2016c.mrt')
    good_names = [n for n, g in galaxies.items() if len(g['R']) >= 5 and np.max(g['Vobs']) > 10]
    print(f"Galaxies after quality cuts: {len(good_names)}")

    models = OrderedDict([
        ('Holographic', (g_holographic, a_u)),
        ('MOND_simple', (g_mond_simple, a0)),
        ('MOND_standard', (g_mond_standard, a0)),
    ])

    scenarios = [
        ("RAW ERRORS (no floor)",                None,  None),
        ("FLOOR: 2 km/s absolute only",          None,  2.0),
        ("FLOOR: 3% of Vobs",                    3.0,   None),
        ("FLOOR: max(2 km/s, 5% Vobs) [paper]",  5.0,   2.0),
        ("FLOOR: 10% of Vobs",                   10.0,  None),
    ]

    all_stats = {}
    for label, fpct, fabs in scenarios:
        results = run_scenario(galaxies, good_names, models, label, fpct, fabs)
        all_stats[label] = print_results(label, results, models)

    # ── SUMMARY TABLE ──
    print(f"\n\n{'='*90}")
    print(f"  ROBUSTNESS SUMMARY: Holographic median χ²/ν and head-to-head wins across error floors")
    print(f"{'='*90}")
    print(f"{'Error floor':<40} {'Holo med':>9} {'MOND_s med':>10} {'Holo wins':>10}")
    print('-'*75)
    for label, _, _ in scenarios:
        s = all_stats.get(label, {})
        h_med = s.get('Holographic', {}).get('med', 0)
        m_med = s.get('MOND_simple', {}).get('med', 0)
        # Re-derive wins from stored data (approximate from medians)
        print(f"{label:<40} {h_med:>9.2f} {m_med:>10.2f}")

    print(f"\nDone.")

if __name__ == '__main__':
    main()
