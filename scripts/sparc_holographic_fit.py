#!/usr/bin/env python3
"""
Full SPARC holographic rotation curve analysis.
Reads MassModels_Lelli2016c.mrt directly, fits all 175 galaxies.
"""
import os, sys, warnings
import numpy as np
from scipy.optimize import minimize_scalar, minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import OrderedDict

warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════
C_MS     = 2.998e8       # m/s
KPC_M    = 3.086e19      # m per kpc
G_SI     = 6.674e-11     # m^3 kg^-1 s^-2
MSUN     = 1.989e30      # kg

def au_from_H0(H0):
    """a_u = c*H0/(2*pi), H0 in km/s/Mpc."""
    H0_si = H0 * 1e3 / 3.086e22
    return C_MS * H0_si / (2 * np.pi)

# ═══════════════════════════════════════════════════════════════
# GRAVITY MODELS
# ═══════════════════════════════════════════════════════════════
def g_holographic(gN, au):
    """g_eff = g_N / (1 - exp(-sqrt(|g_N|/a_u)))"""
    gN = np.asarray(gN, dtype=float)
    x = np.sqrt(np.abs(gN) / au)
    with np.errstate(over='ignore', invalid='ignore'):
        denom = 1.0 - np.exp(-x)
    # Taylor for small x: denom ~ x, so g ~ gN/x = sqrt(|gN|*au)
    safe = x > 1e-8
    out = np.where(safe, gN / np.where(safe, denom, 1.0),
                   np.sign(gN + 1e-30) * np.sqrt(np.abs(gN) * au))
    return out

def g_mond_simple(gN, a0):
    """MOND simple: mu(x)=x/(1+x). Algebraic: g = (gN + sqrt(gN^2 + 4|gN|a0))/2"""
    gN = np.asarray(gN, dtype=float)
    return 0.5 * (gN + np.sqrt(gN**2 + 4*np.abs(gN)*a0))

def g_mond_standard(gN, a0):
    """MOND standard: mu(x)=x/sqrt(1+x^2). Iterative solution."""
    gN = np.asarray(gN, dtype=float)
    gN_abs = np.abs(gN)
    g = np.sqrt(gN_abs**2 + gN_abs * a0)  # initial guess
    for _ in range(30):
        x = g / a0
        mu = x / np.sqrt(1 + x**2)
        f = g * mu - gN_abs
        mu_p = 1.0 / (1 + x**2)**1.5
        df = mu + g * mu_p / a0
        g = np.maximum(g - f / (df + 1e-30), 1e-30)
    return g * np.sign(gN + 1e-30)

# ═══════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════
def load_mass_models(filepath):
    """Parse MassModels_Lelli2016c.mrt into dict of galaxies."""
    galaxies = OrderedDict()
    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip('\r\n')
            # Skip header lines
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

            # Try to parse data line
            parts = line.split()
            if len(parts) < 7:
                continue

            try:
                name = parts[0]
                dist = float(parts[1])
                R = float(parts[2])
                Vobs = float(parts[3])
                eV = float(parts[4])
                Vgas = float(parts[5])
                Vdisk = float(parts[6])
                Vbul = float(parts[7]) if len(parts) > 7 else 0.0
            except (ValueError, IndexError):
                continue

            if name not in galaxies:
                galaxies[name] = {'R': [], 'Vobs': [], 'eV': [], 'Vgas': [],
                                  'Vdisk': [], 'Vbul': [], 'dist': dist}
            galaxies[name]['R'].append(R)
            galaxies[name]['Vobs'].append(Vobs)
            galaxies[name]['eV'].append(eV)
            galaxies[name]['Vgas'].append(Vgas)
            galaxies[name]['Vdisk'].append(Vdisk)
            galaxies[name]['Vbul'].append(Vbul)

    # Convert to numpy arrays
    for name in galaxies:
        for key in ['R', 'Vobs', 'eV', 'Vgas', 'Vdisk', 'Vbul']:
            galaxies[name][key] = np.array(galaxies[name][key])
    return galaxies

# ═══════════════════════════════════════════════════════════════
# FITTING
# ═══════════════════════════════════════════════════════════════
def compute_gN_from_V(R_kpc, Vgas, Vdisk, Vbul, Yd, Yb):
    """Compute Newtonian acceleration from velocity components.
    V^2 components are signed: V can be negative (meaning V^2 contribution is negative).
    g_N = V_tot^2 / r in physical units."""
    # V^2 preserving sign: V*|V|
    V2_tot = (Yd * Vdisk * np.abs(Vdisk) +
              Yb * Vbul * np.abs(Vbul) +
              Vgas * np.abs(Vgas))  # (km/s)^2, signed
    r_m = R_kpc * KPC_M
    gN = V2_tot * 1e6 / r_m  # m/s^2
    return gN

def fit_galaxy(gal, model_func, accel, name=""):
    """Fit one galaxy. Returns result dict."""
    R = gal['R']
    Vobs = gal['Vobs']
    eV = gal['eV']
    Vgas = gal['Vgas']
    Vdisk = gal['Vdisk']
    Vbul = gal['Vbul']

    # Minimum error floor: 2 km/s or 5% of Vobs
    eV_use = np.maximum(eV, np.maximum(2.0, 0.05 * np.abs(Vobs)))

    has_bulge = np.any(np.abs(Vbul) > 0.5)

    def chi2(params):
        if has_bulge:
            Yd, Yb = params[0], params[1]
        else:
            Yd = params[0]
            Yb = 0.0

        gN = compute_gN_from_V(R, Vgas, Vdisk, Vbul, Yd, Yb)
        g_eff = model_func(gN, accel)
        Vpred = np.sqrt(np.abs(g_eff) * R * KPC_M) / 1e3  # km/s
        return np.sum(((Vobs - Vpred) / eV_use)**2)

    # Try multiple starting points
    best = None
    starts = [[0.3], [0.5], [0.8], [1.2]]
    if has_bulge:
        starts = [[0.3, 0.3], [0.5, 0.5], [0.5, 0.8], [0.8, 0.5], [1.0, 1.0]]

    bounds = [(0.01, 10.0)]
    if has_bulge:
        bounds = [(0.01, 10.0), (0.01, 10.0)]

    for x0 in starts:
        try:
            res = minimize(chi2, x0, method='L-BFGS-B', bounds=bounds)
            if best is None or res.fun < best.fun:
                best = res
        except:
            continue

    if best is None:
        return None

    Yd = best.x[0]
    Yb = best.x[1] if has_bulge else 0.0
    n_data = len(Vobs)
    n_params = 2 if has_bulge else 1
    dof = max(n_data - n_params, 1)

    gN = compute_gN_from_V(R, Vgas, Vdisk, Vbul, Yd, Yb)
    g_eff = model_func(gN, accel)
    Vpred = np.sqrt(np.abs(g_eff) * R * KPC_M) / 1e3

    return {
        'name': name,
        'chi2': best.fun,
        'dof': dof,
        'chi2r': best.fun / dof,
        'Yd': Yd, 'Yb': Yb,
        'has_bulge': has_bulge,
        'n_data': n_data,
        'Vpred': Vpred,
        'gN': gN,
        'g_eff': g_eff,
        'R': R, 'Vobs': Vobs, 'eV': eV_use,
        'Vgas': Vgas, 'Vdisk': Vdisk, 'Vbul': Vbul,
    }

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════
def main():
    H0 = 67.4
    a_u = au_from_H0(H0)
    a0  = 1.2e-10

    print("="*70)
    print("  SPARC HOLOGRAPHIC ROTATION CURVE ANALYSIS")
    print(f"  H0 = {H0} km/s/Mpc")
    print(f"  a_u = cH0/(2π) = {a_u:.6e} m/s²")
    print(f"  a0 (MOND)      = {a0:.6e} m/s²")
    print(f"  Ratio a_u/a0   = {a_u/a0:.4f}  ({(1-a_u/a0)*100:.1f}% below)")
    print("="*70)

    # Load data
    datafile = '/home/claude/SPARC/MassModels_Lelli2016c.mrt'
    print(f"\nLoading {datafile}...")
    galaxies = load_mass_models(datafile)
    print(f"Loaded {len(galaxies)} galaxies.")

    # Quality cuts
    good_galaxies = OrderedDict()
    for name, gal in galaxies.items():
        if len(gal['R']) >= 5 and np.max(gal['Vobs']) > 10:
            good_galaxies[name] = gal
    print(f"After quality cuts (≥5 points, Vmax > 10 km/s): {len(good_galaxies)} galaxies.\n")

    # Models
    models = OrderedDict([
        ('Holographic', (g_holographic, a_u)),
        ('MOND_simple', (g_mond_simple, a0)),
        ('MOND_standard', (g_mond_standard, a0)),
    ])

    all_results = {m: [] for m in models}
    outdir = '/home/claude/results'
    os.makedirs(outdir, exist_ok=True)

    # Fit all galaxies
    for i, (gname, gal) in enumerate(good_galaxies.items()):
        for mname, (mfunc, accel) in models.items():
            res = fit_galaxy(gal, mfunc, accel, name=gname)
            if res is not None:
                all_results[mname].append(res)
        if (i+1) % 25 == 0:
            print(f"  Fitted {i+1}/{len(good_galaxies)}...")

    print(f"\n{'Model':<18} {'Med χ²/ν':>9} {'Mean χ²/ν':>10} "
          f"{'σ':>7} {'<1':>6} {'<2':>6} {'<5':>6} {'Med Υd':>7} {'N':>4}")
    print("-"*78)

    stats = {}
    for mname in models:
        rl = all_results[mname]
        if not rl:
            continue
        c2r = np.array([r['chi2r'] for r in rl])
        yds = np.array([r['Yd'] for r in rl])
        med = np.median(c2r)
        mn = np.mean(c2r)
        sd = np.std(c2r)
        f1 = np.mean(c2r < 1)*100
        f2 = np.mean(c2r < 2)*100
        f5 = np.mean(c2r < 5)*100
        my = np.median(yds)
        stats[mname] = {'med': med, 'mean': mn, 'std': sd, 'f1': f1, 'f2': f2, 'f5': f5, 'med_Yd': my}
        print(f"{mname:<18} {med:>9.2f} {mn:>10.2f} {sd:>7.2f} "
              f"{f1:>5.1f}% {f2:>5.1f}% {f5:>5.1f}% {my:>7.3f} {len(rl):>4}")

    # ═══════════════════════════════════════════════════════════
    # SAVE CSV
    # ═══════════════════════════════════════════════════════════
    with open(f'{outdir}/fit_summary.csv', 'w') as f:
        f.write("galaxy,model,chi2,dof,chi2_red,Yd,Yb,n_data,has_bulge\n")
        for mname in models:
            for r in all_results[mname]:
                f.write(f"{r['name']},{mname},{r['chi2']:.4f},{r['dof']},"
                        f"{r['chi2r']:.4f},{r['Yd']:.4f},{r['Yb']:.4f},"
                        f"{r['n_data']},{r['has_bulge']}\n")

    # Save global stats
    with open(f'{outdir}/global_statistics.txt', 'w') as f:
        f.write(f"SPARC Holographic Rotation Curve Analysis\n")
        f.write(f"H0 = {H0} km/s/Mpc\n")
        f.write(f"a_u = cH0/(2π) = {a_u:.6e} m/s²\n")
        f.write(f"a0 (MOND)      = {a0:.6e} m/s²\n")
        f.write(f"Ratio a_u/a0   = {a_u/a0:.4f}\n")
        f.write(f"Galaxies fitted: {len(good_galaxies)}\n\n")
        f.write(f"{'Model':<18} {'Med χ²/ν':>9} {'Mean':>7} {'<1':>6} {'<2':>6} {'<5':>6} {'Med Υd':>7}\n")
        for mname, s in stats.items():
            f.write(f"{mname:<18} {s['med']:>9.2f} {s['mean']:>7.2f} "
                    f"{s['f1']:>5.1f}% {s['f2']:>5.1f}% {s['f5']:>5.1f}% {s['med_Yd']:>7.3f}\n")

    # ═══════════════════════════════════════════════════════════
    # FIGURES
    # ═══════════════════════════════════════════════════════════
    plt.rcParams.update({'font.size': 9, 'font.family': 'serif', 'figure.dpi': 150,
                         'savefig.dpi': 300, 'axes.linewidth': 0.7})

    colors = {'Holographic': '#2563eb', 'MOND_simple': '#dc2626', 'MOND_standard': '#16a34a'}
    labels = {'Holographic': r'Holographic ($a_u = cH_0/2\pi$)',
              'MOND_simple': r'MOND simple ($a_0$)',
              'MOND_standard': r'MOND standard ($a_0$)'}

    # ── Fig 1: Υ_disk histogram ────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5), sharey=True)
    for ax, mname in zip(axes, models):
        yds = [r['Yd'] for r in all_results[mname]]
        ax.hist(yds, bins=np.linspace(0, 3, 31), color=colors[mname],
                alpha=0.7, edgecolor='white', lw=0.5)
        m = np.median(yds)
        ax.axvline(m, color='k', ls='--', lw=1, label=f'Median={m:.2f}')
        ax.axvspan(0.2, 0.8, color='gray', alpha=0.12, label='SPS range')
        ax.set_xlabel(r'$\Upsilon_{\rm disk}$ [$M_\odot/L_\odot$]')
        ax.set_title(labels[mname], fontsize=8)
        ax.legend(fontsize=7)
        ax.set_xlim(0, 3)
    axes[0].set_ylabel('Count')
    fig.tight_layout()
    fig.savefig(f'{outdir}/fig1_upsilon_hist.pdf', bbox_inches='tight')
    fig.savefig(f'{outdir}/fig1_upsilon_hist.png', bbox_inches='tight')
    plt.close(fig)
    print("\n  Fig 1 saved.")

    # ── Fig 2: χ² distribution ─────────────────────────────────
    fig, ax = plt.subplots(figsize=(6, 4))
    for mname in models:
        c2r = sorted([r['chi2r'] for r in all_results[mname]])
        cdf = np.arange(1, len(c2r)+1) / len(c2r)
        ax.plot(c2r, cdf, color=colors[mname], lw=1.5, label=labels[mname])
    ax.set_xlabel(r'Reduced $\chi^2$')
    ax.set_ylabel('Cumulative fraction')
    ax.set_xlim(0, 15)
    ax.axvline(1, color='gray', ls=':', lw=0.5)
    ax.axvline(2, color='gray', ls=':', lw=0.5)
    ax.legend(fontsize=8)
    ax.set_title(r'$\chi^2_\nu$ Distribution — 175 SPARC Galaxies', fontsize=10)
    fig.tight_layout()
    fig.savefig(f'{outdir}/fig2_chi2_cdf.pdf', bbox_inches='tight')
    fig.savefig(f'{outdir}/fig2_chi2_cdf.png', bbox_inches='tight')
    plt.close(fig)
    print("  Fig 2 saved.")

    # ── Fig 3: Radial Acceleration Relation ────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, mname, accel in zip(axes, ['Holographic', 'MOND_simple'], [a_u, a0]):
        gbar_all, gobs_all = [], []
        for r in all_results[mname]:
            gN = r['gN']
            gobs = (r['Vobs']*1e3)**2 / (r['R']*KPC_M)
            ok = (gN > 0) & (gobs > 0)
            gbar_all.extend(gN[ok])
            gobs_all.extend(gobs[ok])

        gbar_all = np.array(gbar_all)
        gobs_all = np.array(gobs_all)

        ax.scatter(np.log10(gbar_all), np.log10(gobs_all),
                   s=0.2, alpha=0.1, color=colors[mname], rasterized=True)

        # 1:1 line
        xr = np.linspace(-13, -8, 200)
        ax.plot(xr, xr, 'k-', lw=0.6, label='1:1')

        # Model line
        gt = 10**xr
        if mname == 'Holographic':
            gp = g_holographic(gt, accel)
        else:
            gp = g_mond_simple(gt, accel)
        ax.plot(xr, np.log10(gp), '--', color=colors[mname], lw=1.5, label=labels[mname])

        ax.set_xlabel(r'$\log_{10}\, g_{\rm bar}$ [m/s$^2$]')
        ax.set_ylabel(r'$\log_{10}\, g_{\rm obs}$ [m/s$^2$]')
        ax.set_xlim(-12.5, -8.5)
        ax.set_ylim(-12.5, -8.5)
        ax.set_aspect('equal')
        ax.legend(fontsize=7)
        ax.set_title(labels[mname], fontsize=8)

    fig.tight_layout()
    fig.savefig(f'{outdir}/fig3_rar.pdf', bbox_inches='tight')
    fig.savefig(f'{outdir}/fig3_rar.png', bbox_inches='tight')
    plt.close(fig)
    print("  Fig 3 saved.")

    # ── Fig 4: Baryonic Tully-Fisher ───────────────────────────
    fig, ax = plt.subplots(figsize=(5.5, 5))
    for mname, marker, ms in [('Holographic', 'o', 2.5), ('MOND_simple', 's', 2)]:
        accel = a_u if mname == 'Holographic' else a0
        logMb_list, logvf_list = [], []
        for r in all_results[mname]:
            vf = np.median(r['Vobs'][-max(3, len(r['Vobs'])//5):])
            if vf < 15:
                continue
            # Mb from BTFR: v^4 = G*Mb*a => Mb = v^4/(G*a)
            Mb = (vf*1e3)**4 / (G_SI * accel)
            logvf_list.append(np.log10(vf))
            logMb_list.append(np.log10(Mb / MSUN))

        ax.scatter(logvf_list, logMb_list, s=ms**2, alpha=0.4,
                   color=colors[mname], marker=marker, label=labels[mname])

    # Predicted line: log Mb = 4 log v - log(G*a_u) with v in km/s
    vr = np.linspace(1.0, 2.7, 100)
    logMb_pred = 4*vr + np.log10((1e3)**4 / (G_SI * a_u * MSUN))
    ax.plot(vr, logMb_pred, 'b-', lw=1.5, label=r'$v^4 = GM_b a_u$ (zero params)')

    ax.set_xlabel(r'$\log_{10}\, v_f$ [km/s]')
    ax.set_ylabel(r'$\log_{10}\, M_b / M_\odot$')
    ax.legend(fontsize=7)
    ax.set_title('Baryonic Tully–Fisher Relation', fontsize=10)
    fig.tight_layout()
    fig.savefig(f'{outdir}/fig4_btfr.pdf', bbox_inches='tight')
    fig.savefig(f'{outdir}/fig4_btfr.png', bbox_inches='tight')
    plt.close(fig)
    print("  Fig 4 saved.")

    # ── Fig 5: Galaxy gallery (12 representative) ──────────────
    hres = all_results['Holographic']
    sorted_h = sorted(hres, key=lambda r: r['chi2r'])
    n_total = len(sorted_h)
    # Pick 12: best 4, median 4, worst 4
    picks = []
    for frac in [0.02, 0.08, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.92, 0.98]:
        picks.append(sorted_h[min(int(frac * n_total), n_total-1)])

    fig, axes = plt.subplots(3, 4, figsize=(13, 9))
    for idx, ax in enumerate(axes.flat):
        if idx >= len(picks):
            ax.set_visible(False)
            continue
        r = picks[idx]
        ax.errorbar(r['R'], r['Vobs'], yerr=r['eV'], fmt='o', ms=2, color='k',
                     elinewidth=0.4, capsize=0, label='Data', zorder=3)
        ax.plot(r['R'], r['Vpred'], '-', color='#2563eb', lw=1.5, label='Holographic')

        # Baryonic only
        Vbar = np.sqrt(np.abs(
            r['Yd']*r['Vdisk']*np.abs(r['Vdisk']) +
            r['Yb']*r['Vbul']*np.abs(r['Vbul']) +
            r['Vgas']*np.abs(r['Vgas'])
        ))
        ax.plot(r['R'], Vbar, ':', color='gray', lw=0.8, label='Baryons')

        ax.set_title(f"{r['name']}  χ²/ν={r['chi2r']:.1f}  Υ={r['Yd']:.2f}", fontsize=7)
        ax.tick_params(labelsize=6)
        if idx >= 8:
            ax.set_xlabel('R [kpc]', fontsize=7)
        if idx % 4 == 0:
            ax.set_ylabel('V [km/s]', fontsize=7)
        if idx == 0:
            ax.legend(fontsize=5)

    fig.suptitle('Holographic Fits: Best → Worst (left to right, top to bottom)', fontsize=10, y=1.01)
    fig.tight_layout()
    fig.savefig(f'{outdir}/fig5_gallery.pdf', bbox_inches='tight')
    fig.savefig(f'{outdir}/fig5_gallery.png', bbox_inches='tight')
    plt.close(fig)
    print("  Fig 5 saved.")

    # ── Fig 6: Interpolation function comparison ───────────────
    fig, ax = plt.subplots(figsize=(6, 4))
    x = np.linspace(-2, 3, 500)
    gt = 10**x

    # mu = g_N / g_eff for each model
    mu_h = (gt * a_u) / g_holographic(gt * a_u, a_u)
    mu_s = (gt * a0) / g_mond_simple(gt * a0, a0)
    mu_st = (gt * a0) / g_mond_standard(gt * a0, a0)

    ax.plot(x, mu_h, '-', color='#2563eb', lw=2, label='Holographic')
    ax.plot(x, mu_s, '--', color='#dc2626', lw=1.5, label='MOND simple')
    ax.plot(x, mu_st, '-.', color='#16a34a', lw=1.5, label='MOND standard')
    ax.axhline(1, color='gray', ls=':', lw=0.5)
    ax.set_xlabel(r'$\log_{10}(g_N / a)$')
    ax.set_ylabel(r'$\mu(g_N/a)$')
    ax.legend(fontsize=9)
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlim(-2, 3)
    ax.set_title('Interpolation Function Comparison', fontsize=10)
    fig.tight_layout()
    fig.savefig(f'{outdir}/fig6_transition.pdf', bbox_inches='tight')
    fig.savefig(f'{outdir}/fig6_transition.png', bbox_inches='tight')
    plt.close(fig)
    print("  Fig 6 saved.")

    # ── Fig 7: Δχ² holographic vs MOND ─────────────────────────
    fig, ax = plt.subplots(figsize=(6, 4))
    # Match galaxies between holographic and MOND simple
    h_dict = {r['name']: r['chi2r'] for r in all_results['Holographic']}
    m_dict = {r['name']: r['chi2r'] for r in all_results['MOND_simple']}
    common = sorted(set(h_dict.keys()) & set(m_dict.keys()))

    h_vals = np.array([h_dict[n] for n in common])
    m_vals = np.array([m_dict[n] for n in common])
    delta = h_vals - m_vals  # positive = holographic worse

    ax.hist(delta, bins=np.linspace(-5, 5, 51), color='#7c3aed', alpha=0.7, edgecolor='white', lw=0.5)
    ax.axvline(0, color='k', ls='-', lw=1)
    ax.axvline(np.median(delta), color='red', ls='--', lw=1.2,
               label=f'Median Δ = {np.median(delta):.2f}')
    n_better = np.sum(delta < 0)
    n_worse = np.sum(delta > 0)
    ax.set_xlabel(r'$\chi^2_\nu$(Holographic) $-$ $\chi^2_\nu$(MOND simple)')
    ax.set_ylabel('Count')
    ax.set_title(f'Holographic better: {n_better}/{len(common)}, '
                 f'MOND better: {n_worse}/{len(common)}', fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f'{outdir}/fig7_delta_chi2.pdf', bbox_inches='tight')
    fig.savefig(f'{outdir}/fig7_delta_chi2.png', bbox_inches='tight')
    plt.close(fig)
    print("  Fig 7 saved.")

    print(f"\n{'='*70}")
    print(f"  ALL DONE. Results in {outdir}/")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
