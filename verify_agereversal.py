# -*- coding: utf-8 -*-
"""
裁定实验：改进方案(年龄反转)的数字地基核查
==================================================
目标: Stage14(del代码, cogtot) 报 younger β=-0.261(p=.009) / older β=+0.107(p=.344) / Wald p=.015
      vs 09C(cog27) 报 older unweighted β=+0.344(p=.004) 显著
      两者结局列不同(cogtot vs cog27)。用当前 Stage2 parquet 同口径复跑，
      裁定 age-reversal 在稿件主分析结局(cog27)上是否成立，以及差异来源。
口径: 与 Stage14 完全一致: 每人 wave3-15 首次观测的 agey_b 为 baseline age,
      中位数(59)切 younger/older; TWFE 双向demean + hhidpn cluster; x=[copd_ever]+COVS
COVS (Stage14): agey_b,female,raedyrs,race_white,race_black,hispanic,married,
                smoken,cesd,bmi,hibpe,diabe,stroke
同时用 cog27 与 cogtot 各跑一遍,对比。
"""
import numpy as np, pandas as pd, time
from scipy import stats

PANEL = r"D:\HRSCOPD\数据分析\COPD与认知\投稿代码\code\Pipeline_Config_Refactor\输出结果\hrs_panel_stage2.parquet"

def twoway_demean(df, y_col, x_cols, entity="hhidpn", time="wave", max_iter=300, tol=1e-9):
    cols = [entity, time, y_col] + x_cols
    d = df[[c for c in cols if c in df.columns]].copy()
    for c in [y_col] + x_cols:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna()
    if len(d) < 50: return None
    eid, tid = d[entity].values, d[time].values
    Y = d[y_col].values.astype(float)
    X = d[x_cols].values.astype(float)
    e_codes, e_levels = pd.factorize(eid)
    t_codes, t_levels = pd.factorize(tid)
    n_e, n_t = len(e_levels), len(t_levels)
    def _demean_col(v):
        v = v - v.mean()
        for _ in range(max_iter):
            v0 = v.copy()
            e_sum = np.bincount(e_codes, weights=v, minlength=n_e)
            e_cnt = np.bincount(e_codes, minlength=n_e)
            v = v - (e_sum / np.maximum(e_cnt, 1))[e_codes]
            t_sum = np.bincount(t_codes, weights=v, minlength=n_t)
            t_cnt = np.bincount(t_codes, minlength=n_t)
            v = v - (t_sum / np.maximum(t_cnt, 1))[t_codes]
            if np.max(np.abs(v - v0)) < tol: break
        return v
    Yd = _demean_col(Y.copy())
    Xd = np.column_stack([_demean_col(X[:, j].copy()) for j in range(X.shape[1])])
    return Yd, Xd, eid, len(d)

def clustered_se(Xd, resid, clusters):
    n, k = Xd.shape
    XtXi = np.linalg.pinv(Xd.T @ Xd)
    G = len(np.unique(clusters))
    meat = np.zeros((k, k))
    for cl in np.unique(clusters):
        m = clusters == cl
        s = Xd[m].T @ resid[m]
        meat += np.outer(s, s)
    corr = (G / (G - 1)) * ((n - 1) / (n - k))
    return np.sqrt(np.diag(corr * XtXi @ meat @ XtXi))

def run_fe(df, outcome, treatment, covars=None, entity="hhidpn", time="wave"):
    covars = [v for v in (covars or []) if v in df.columns and v != treatment]
    xcols = [treatment] + covars
    res = twoway_demean(df, outcome, xcols, entity, time)
    if res is None: return None
    Yd, Xd, eid, n = res
    k = Xd.shape[1]
    if n < 100 or n <= k: return None
    betas = np.linalg.lstsq(Xd, Yd, rcond=None)[0]
    resid = Yd - Xd @ betas
    n_ent = len(np.unique(eid))
    try:
        ses = clustered_se(Xd, resid, eid)
    except Exception:
        sigma2 = (resid**2).sum() / max(n - k, 1)
        ses = np.sqrt(np.diag(sigma2 * np.linalg.pinv(Xd.T @ Xd)))
    coef, se = betas[0], ses[0]
    tval = coef / se if se > 0 else np.nan
    df_ = max(n - n_ent - k, 1)
    pval = 2 * (1 - stats.t.cdf(abs(tval), df=df_))
    return {"coef": coef, "se": se, "p": pval, "n_obs": n, "n_ent": n_ent}

def wald_compare(r1, r2):
    d = r1["coef"] - r2["coef"]
    se = np.sqrt(r1["se"]**2 + r2["se"]**2)
    z = d / se
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return d, z, p

t0 = time.time()
print("Load panel ...")
panel = pd.read_parquet(PANEL)
print(f"  panel {panel.shape}, individuals {panel['hhidpn'].nunique():,}, waves {panel['wave'].min()}-{panel['wave'].max()}")
print(f"  cog27 nonnull={panel['cog27'].notna().sum():,}  cogtot nonnull={panel['cogtot'].notna().sum():,}")

COVS = ["agey_b","female","raedyrs","race_white","race_black","hispanic","married",
        "smoken","cesd","bmi","hibpe","diabe","stroke"]
COVS = [c for c in COVS if c in panel.columns]

# baseline age = 每人 wave3-15 首次观测 agey_b
sub = panel[panel["wave"].between(3, 15)].sort_values(["hhidpn", "wave"])
bl = sub.groupby("hhidpn")["agey_b"].first().rename("baseline_age")
panel = panel.merge(bl, on="hhidpn", how="left")
med = panel["baseline_age"].median()
print(f"  baseline-age median = {med:.1f}")
panel["age_group"] = np.where(panel["baseline_age"] < med, "younger", "older")
print("  n per group (person-wave in 3-15):")
print(panel[panel["wave"].between(3,15)].groupby("age_group")["hhidpn"].nunique())

for cog_col in ["cogtot", "cog27"]:
    print(f"\n===== outcome = {cog_col} =====")
    rows = {}
    for grp in ["younger", "older"]:
        d = panel[(panel["age_group"] == grp) & (panel["wave"].between(3, 15))]
        r = run_fe(d, cog_col, "copd_ever", COVS)
        rows[grp] = r
        if r:
            stars = "***" if r["p"]<.001 else "**" if r["p"]<.01 else "*" if r["p"]<.05 else "n.s."
            print(f"  {grp:8s}  beta={r['coef']:+.4f}  SE={r['se']:.4f}  p={r['p']:.4f} {stars}  n={r['n_obs']:,}  indiv={r['n_ent']:,}")
    if rows["younger"] and rows["older"]:
        d, z, p = wald_compare(rows["younger"], rows["older"])
        print(f"  Wald younger-vs-older: Delta={d:+.4f}  z={z:+.3f}  p={p:.4f}")
print(f"\nDone in {time.time()-t0:.1f}s")
