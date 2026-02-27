#!/usr/bin/env python3
"""
Isolation test: Is it the acceleration scale or the interpolation function?
2x2: {Holo, MOND_simple} x {a_u, a_0}
"""
import os, warnings
import numpy as np
from scipy.optimize import minimize
from collections import OrderedDict
warnings.filterwarnings('ignore')

C_MS = 2.998e8; KPC_M = 3.086e19

def au_from_H0(H0):
    return C_MS * H0 * 1e3 / 3.086e22 / (2*np.pi)

def g_holographic(gN, a):
    gN = np.asarray(gN, dtype=float)
    x = np.sqrt(np.abs(gN)/a)
    with np.errstate(over='ignore', invalid='ignore'):
        d = 1.0 - np.exp(-x)
    safe = x > 1e-8
    return np.where(safe, gN/np.where(safe,d,1.0), np.sign(gN+1e-30)*np.sqrt(np.abs(gN)*a))

def g_mond_simple(gN, a):
    gN = np.asarray(gN, dtype=float)
    return 0.5*(gN + np.sqrt(gN**2 + 4*np.abs(gN)*a))

def load_data(fp):
    gals = OrderedDict()
    with open(fp) as f:
        for line in f:
            line = line.rstrip('\r\n')
            if any(line.startswith(x) for x in ['Title','Author','Table','=','Byte','-','Note']) or \
               line.strip()=='' or any(x in line for x in ['Label','Explanations','Units','uncertainty','factor','solMass','and Accurate']) or \
               (line.startswith(' ') and ('Bytes' in line or 'Format' in line)):
                continue
            p = line.split()
            if len(p)<7: continue
            try:
                nm=p[0]; R=float(p[2]); Vo=float(p[3]); eV=float(p[4])
                Vg=float(p[5]); Vd=float(p[6]); Vb=float(p[7]) if len(p)>7 else 0.0
            except: continue
            if nm not in gals:
                gals[nm]={'R':[],'Vobs':[],'eV':[],'Vgas':[],'Vdisk':[],'Vbul':[]}
            gals[nm]['R'].append(R); gals[nm]['Vobs'].append(Vo); gals[nm]['eV'].append(eV)
            gals[nm]['Vgas'].append(Vg); gals[nm]['Vdisk'].append(Vd); gals[nm]['Vbul'].append(Vb)
    for n in gals:
        for k in gals[n]: gals[n][k]=np.array(gals[n][k])
    return gals

def fit_one(gal, mfunc, accel):
    R,Vo,eV=gal['R'],gal['Vobs'],gal['eV']
    Vg,Vd,Vb=gal['Vgas'],gal['Vdisk'],gal['Vbul']
    eV_use=np.maximum(eV, np.maximum(2.0, 0.05*np.abs(Vo)))
    hb=np.any(np.abs(Vb)>0.5)
    def chi2(par):
        Yd=par[0]; Yb=par[1] if hb else 0.0
        if Yd<0.01 or Yd>10: return 1e10
        if hb and (Yb<0.01 or Yb>10): return 1e10
        V2=Yd*Vd*np.abs(Vd)+Yb*Vb*np.abs(Vb)+Vg*np.abs(Vg)
        gN=V2*1e6/(R*KPC_M)
        ge=mfunc(gN,accel)
        Vp=np.sqrt(np.abs(ge)*R*KPC_M)/1e3
        return np.sum(((Vo-Vp)/eV_use)**2)
    sts=[[0.3],[0.5],[0.8],[1.2]]
    bds=[(0.01,10.0)]
    if hb: sts=[[.3,.3],[.5,.5],[.5,.8],[.8,.5],[1,1]]; bds=[(0.01,10),(0.01,10)]
    best=None
    for x0 in sts:
        try:
            r=minimize(chi2,x0,method='L-BFGS-B',bounds=bds)
            if best is None or r.fun<best.fun: best=r
        except: pass
    if best is None: return None
    dof=max(len(Vo)-(2 if hb else 1),1)
    return {'chi2r':best.fun/dof,'Yd':best.x[0]}

def h2h(rA, rB):
    common=sorted(set(rA.keys())&set(rB.keys()))
    if not common: return 0,0,0,0
    dA=np.array([rA[n]['chi2r'] for n in common])
    dB=np.array([rB[n]['chi2r'] for n in common])
    d=dA-dB
    return len(common), int(np.sum(d<0)), int(np.sum(d>0)), float(np.median(d))

def main():
    a_u=au_from_H0(67.4); a_0=1.2e-10
    print(f"a_u = {a_u:.4e}  |  a_0 = {a_0:.4e}  |  ratio = {a_u/a_0:.4f}\n")

    gals=load_data('/home/claude/SPARC/MassModels_Lelli2016c.mrt')
    good=[n for n,g in gals.items() if len(g['R'])>=5 and np.max(g['Vobs'])>10]
    print(f"Galaxies: {len(good)}\n")

    cfgs=OrderedDict([
        ('Holo+a_u',   (g_holographic, a_u)),
        ('Holo+a_0',   (g_holographic, a_0)),
        ('MONDs+a_0',  (g_mond_simple, a_0)),
        ('MONDs+a_u',  (g_mond_simple, a_u)),
    ])

    res={k:{} for k in cfgs}
    for i,gn in enumerate(good):
        for cn,(fn,ac) in cfgs.items():
            r=fit_one(gals[gn],fn,ac)
            if r: res[cn][gn]=r
        if (i+1)%50==0: print(f"  {i+1}/{len(good)}...")

    print(f"\n{'Config':<16} {'Med χ²/ν':>9} {'Mean':>7} {'<1':>6} {'<2':>6} {'Med Υd':>7}")
    print("-"*52)
    for cn in cfgs:
        c2=np.array([v['chi2r'] for v in res[cn].values()])
        yd=np.array([v['Yd'] for v in res[cn].values()])
        print(f"{cn:<16} {np.median(c2):>9.3f} {np.mean(c2):>7.2f} "
              f"{np.mean(c2<1)*100:>5.1f}% {np.mean(c2<2)*100:>5.1f}% {np.median(yd):>7.3f}")

    print(f"\n{'='*70}")
    print(f"  HEAD-TO-HEAD ISOLATION")
    print(f"{'='*70}")

    # Scale effect
    n,wA,wB,md = h2h(res['Holo+a_u'], res['Holo+a_0'])
    print(f"\n  Q1: SAME FUNCTION (Eq.2) — a_u vs a_0")
    print(f"      a_u wins: {wA}/{n} ({wA/n*100:.1f}%)   a_0 wins: {wB}/{n} ({wB/n*100:.1f}%)   med Δ={md:+.4f}")

    n,wA,wB,md = h2h(res['MONDs+a_u'], res['MONDs+a_0'])
    print(f"\n  Q2: SAME FUNCTION (MOND simple) — a_u vs a_0")
    print(f"      a_u wins: {wA}/{n} ({wA/n*100:.1f}%)   a_0 wins: {wB}/{n} ({wB/n*100:.1f}%)   med Δ={md:+.4f}")

    # Function effect
    n,wA,wB,md = h2h(res['Holo+a_u'], res['MONDs+a_u'])
    print(f"\n  Q3: SAME SCALE (a_u) — Eq.2 vs MOND simple")
    print(f"      Eq.2 wins: {wA}/{n} ({wA/n*100:.1f}%)   MOND wins: {wB}/{n} ({wB/n*100:.1f}%)   med Δ={md:+.4f}")

    n,wA,wB,md = h2h(res['Holo+a_0'], res['MONDs+a_0'])
    print(f"\n  Q4: SAME SCALE (a_0) — Eq.2 vs MOND simple")
    print(f"      Eq.2 wins: {wA}/{n} ({wA/n*100:.1f}%)   MOND wins: {wB}/{n} ({wB/n*100:.1f}%)   med Δ={md:+.4f}")

    # Overall
    n,wA,wB,md = h2h(res['Holo+a_u'], res['MONDs+a_0'])
    print(f"\n  PAPER COMPARISON: Holo+a_u vs MOND_s+a_0")
    print(f"      Holo+a_u wins: {wA}/{n} ({wA/n*100:.1f}%)   med Δ={md:+.4f}")

    # Summary
    print(f"\n{'='*70}")
    print(f"  VERDICT")
    print(f"{'='*70}")

    n1,s1,_,_ = h2h(res['Holo+a_u'], res['Holo+a_0'])
    n2,s2,_,_ = h2h(res['MONDs+a_u'], res['MONDs+a_0'])
    n3,f1,_,_ = h2h(res['Holo+a_u'], res['MONDs+a_u'])
    n4,f2,_,_ = h2h(res['Holo+a_0'], res['MONDs+a_0'])

    scale_pct = ((s1/n1 + s2/n2)/2)*100
    func_pct = ((f1/n3 + f2/n4)/2)*100
    print(f"\n  Average scale effect (a_u vs a_0):     {scale_pct:.1f}% win rate for a_u")
    print(f"  Average function effect (Eq.2 vs MOND): {func_pct:.1f}% win rate for Eq.2")

    if scale_pct > func_pct:
        print(f"\n  >> THE SCALE IS DOING MORE WORK THAN THE FUNCTION. <<")
    elif func_pct > scale_pct:
        print(f"\n  >> THE FUNCTION IS DOING MORE WORK THAN THE SCALE. <<")
    else:
        print(f"\n  >> BOTH CONTRIBUTE EQUALLY. <<")

if __name__=='__main__':
    main()
