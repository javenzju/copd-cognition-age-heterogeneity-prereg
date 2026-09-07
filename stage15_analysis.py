# -*- coding: utf-8 -*-
"""
方案A / Stage 15 综合分析（年龄反转方案的可执行剩余工作）
================================================================
用户已确认：JGSA 无在审稿(排除解除)；OSF 预注册 + medRxiv/bioRxiv/AAIC 检索 + 预印本 暂缓。
本脚本完成其余全部可分析工作：

  Stage 15A : COPD x age-group x DunedinPACE 三向交互 TWFE（仅 pace 支线，telomere 支线放弃）
              + 最小可检测 Δβ（power 80%）作为功效陈述
  Stage 15B : 分年龄层 joint biomarker mediation (NfL + pTau181 + Homocysteine)，方向开放
              （不再预设"年轻组更强"；直接报 two strata 的真实结果）
  Telomere  : 全样本 telomere 效应修饰（cog27，不做年龄分层），复跑整合既有 Wald p=.041

所有模型统一使用稿件主结局 cog27 + FULL 协变量规格（含 agey_b），与 Step1 决策一致。
复用 Stage10B / Stage9 的辅助函数（twoway_demean / clustered_se / 并行多重中介 bootstrap）。
"""
import warnings, numpy as np, pandas as pd
import pyreadstat
from pathlib import Path
from scipy import stats
warnings.filterwarnings("ignore")

PANEL = Path(r"D:\HRSCOPD\数据分析\COPD与认知\投稿代码\code\Pipeline_Config_Refactor\输出结果\hrs_panel_stage2.parquet")
POA   = Path(r"D:\HRS\PaceOfAging\HRSPoA_Shared.dta")
HCAP  = Path(r"D:\HRS\HCAPbioPilot\HCAPBIOPILOTA_R.dta")
VBS   = Path(r"D:\HRS\HRS2016VBSsubs\vbs16ss.dta")
OUT   = Path(r"D:\HRSCOPD\数据分析\COPD与认知\方案A")
OUT.mkdir(parents=True, exist_ok=True)

Z975, Z80 = 1.959964, 0.841621
# ───────────────────────────── 辅助函数 ─────────────────────────────
def load_dta(p, cols=None):
    if cols:
        _, meta = pyreadstat.read_dta(str(p), metadataonly=True)
        avail = {c.lower() for c in meta.column_names}
        use = [c for c in cols if c.lower() in avail]
        df, _ = pyreadstat.read_dta(str(p), usecols=use)
    else:
        df, _ = pyreadstat.read_dta(str(p))
    df.columns = df.columns.str.lower().str.strip()
    return df

def make_hhidpn(df):
    if "hhidpn" in df.columns:
        df["hhidpn"] = pd.to_numeric(df["hhidpn"], errors="coerce").astype("Int64")
        return df
    hh = next((c for c in df.columns if c == "hhid"), None)
    pn = next((c for c in df.columns if c == "pn"), None)
    h = df[hh].astype(str).str.strip().str.zfill(6)
    p = df[pn].astype(str).str.strip().str.zfill(3)
    df["hhidpn"] = (h + p).astype(np.int64)
    return df

def twoway_demean(y, X, eid, tid, max_iter=300, tol=1e-9):
    e_codes, e_levels = pd.factorize(eid); t_codes, t_levels = pd.factorize(tid)
    n_e, n_t = len(e_levels), len(t_levels)
    def _dm(v):
        v = v - v.mean()
        for _ in range(max_iter):
            v0 = v.copy()
            v = v - (np.bincount(e_codes, weights=v, minlength=n_e)/np.bincount(e_codes, minlength=n_e))[e_codes]
            v = v - (np.bincount(t_codes, weights=v, minlength=n_t)/np.bincount(t_codes, minlength=n_t))[t_codes]
            if np.max(np.abs(v - v0)) < tol: break
        return v
    Yd = _dm(y.astype(float).copy())
    Xd = np.column_stack([_dm(X[:, j].copy()) for j in range(X.shape[1])])
    return Yd, Xd

def clustered_se(Xd, resid, eid):
    n, k = Xd.shape
    XtXi = np.linalg.pinv(Xd.T @ Xd)
    G = len(np.unique(eid))
    meat = np.zeros((k, k))
    for cl in np.unique(eid):
        m = eid == cl
        s = Xd[m].T @ resid[m]
        meat += np.outer(s, s)
    corr = (G/(G-1)) * ((n-1)/(n-k))
    return np.sqrt(np.diag(corr * XtXi @ meat @ XtXi))

def fe_terms(df, ycol, term_cols, covars, entity="hhidpn", time="wave"):
    """TWFE 后 OLS + 实体聚类 SE；term_cols 为已构造的交互项列名。"""
    need = [ycol, entity, time] + term_cols + [c for c in covars if c in df.columns]
    d = df[[c for c in need if c in df.columns]].copy()
    for c in [ycol]+term_cols:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna()
    if len(d) < 200: return None
    Y = d[ycol].values.astype(float)
    X = d[term_cols + [c for c in covars if c in d.columns]].values.astype(float)
    eid = d[entity].values; tid = d[time].values
    Yd, Xd = twoway_demean(Y, X, eid, tid)
    n, k = Xd.shape
    if n <= k: return None
    betas = np.linalg.lstsq(Xd, Yd, rcond=None)[0]
    resid = Yd - Xd @ betas
    n_ent = len(np.unique(eid))
    try:
        ses = clustered_se(Xd, resid, eid)
    except Exception:
        ses = np.sqrt(np.diag((resid@resid)/(n-k) * np.linalg.pinv(Xd.T@Xd)))
    out = {}
    for i, name in enumerate(term_cols + [c for c in covars if c in d.columns]):
        coef, se = betas[i], ses[i]
        t = coef/se if se>0 else np.nan
        dfree = max(n - n_ent - k, 1)
        p = 2*(1-stats.t.cdf(abs(t), df=dfree))
        out[name] = dict(coef=coef, se=se, p=p,
                         ci_lo=coef-1.96*se, ci_hi=coef+1.96*se, n_obs=n, n_ent=n_ent)
    return out

def joint_mediation(d, meds, covs, n_boot=2000, seed=20260907):
    """单次并行多重中介 + bootstrap；d 必须已 dropna 完整。返回 dict。"""
    cols = [COG, COPD] + meds + [c for c in covs if c in d.columns]
    dd = d[cols].apply(pd.to_numeric, errors="coerce").dropna()
    n = len(dd)
    if n < 30: return None
    x_c = [c for c in covs if c in dd.columns]
    def _est(dd):
        a = {}
        for m in meds:
            X = np.column_stack([np.ones(n), dd[COPD].values] + [dd[c].values for c in x_c])
            b,_,_,_ = np.linalg.lstsq(X, dd[m].values, rcond=None); a[m]=b[1]
        X2 = np.column_stack([np.ones(n), dd[COPD].values] + [dd[m].values for m in meds] + [dd[c].values for c in x_c])
        b2,_,_,_ = np.linalg.lstsq(X2, dd[COG].values, rcond=None)
        c_prime = b2[1]; b_paths = {m: b2[2+i] for i,m in enumerate(meds)}
        ind = {m: a[m]*b_paths[m] for m in meds}
        X3 = np.column_stack([np.ones(n), dd[COPD].values] + [dd[c].values for c in x_c])
        te,_,_,_ = np.linalg.lstsq(X3, dd[COG].values, rcond=None)
        return dict(joint=sum(ind.values()), partial=ind, total=te[1], direct=c_prime)
    pt = _est(dd)
    rng = np.random.default_rng(seed)
    bj, bp = [], {m:[] for m in meds}
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        r = _est(dd.iloc[idx])
        if r is None: continue
        bj.append(r["joint"]); 
        for m in meds: bp[m].append(r["partial"][m])
    bj = np.array(bj)
    j_lo, j_hi = np.percentile(bj, [2.5,97.5]); j_p = 2*min((bj>0).mean(),(bj<0).mean())
    res = dict(n=n, total=pt["total"], direct=pt["direct"], joint=pt["joint"],
              joint_ci=(j_lo,j_hi), joint_p=j_p, partial={}, partial_p={})
    for m in meds:
        arr = np.array(bp[m]); lo,hi = np.percentile(arr,[2.5,97.5])
        res["partial"][m] = (pt["partial"][m], lo, hi, 2*min((arr>0).mean(),(arr<0).mean()))
    return res

# ════════════════════════════════════════════════════════════════════════
print("LOAD panel + pace + baseline age")
panel = pd.read_parquet(PANEL)
COG, COPD = "cog27", "copd_ever"
COVS = [c for c in ["agey_b","female","raedyrs","race_white","race_black","hispanic",
                    "married","smoken","cesd","bmi","hibpe","diabe","stroke"]
        if c in panel.columns]
poa = make_hhidpn(load_dta(POA, ["hhidpn","paceofaging"]))
poa["pace"] = pd.to_numeric(poa["paceofaging"], errors="coerce")
poa = poa[["hhidpn","pace"]].dropna().drop_duplicates("hhidpn")

# baseline (first-wave) age per person
sub = panel[panel["wave"].between(3,13)].sort_values(["hhidpn","wave"])
base = sub.groupby("hhidpn")["agey_b"].first().rename("baseline_age")
panel = panel.merge(base, left_on="hhidpn", right_index=True, how="left")
panel = panel.merge(poa, on="hhidpn", how="left")
panel["age_grp59"] = (panel["baseline_age"] < 59).astype(float)      # 1=younger
panel["age_grp61"] = (panel["baseline_age"] < panel["baseline_age"].median()).astype(float)
print(f"  panel {panel.shape}, pace merged n_indiv={panel['pace'].notna().sum():,}, "
      f"baseline_age median={panel['baseline_age'].median():.1f}")

# ════════════════════════════════════════════════════════════════════════
# Stage 15A : COPD x age x pace 三向交互 TWFE（仅 pace 支线）
# ════════════════════════════════════════════════════════════════════════
print("\n=== Stage 15A: COPD x age x pace three-way TWFE ===")
rows15a = []
for cut, grpcol in [("age<59", "age_grp59"), ("median-split", "age_grp61")]:
    d = panel.dropna(subset=[COG, COPD, "pace", grpcol] + COVS).copy()
    d["copd"] = d[COPD]
    d["c_ag"]  = d["copd"] * d[grpcol]
    d["c_pa"]  = d["copd"] * d["pace"]
    d["c_agp"] = d["copd"] * d[grpcol] * d["pace"]
    res = fe_terms(d, COG, ["copd","c_ag","c_pa","c_agp"], COVS)
    if res is None:
        print(f"  [{cut}] FAILED (n too small)"); continue
    r3 = res["c_agp"]
    md_beta = (Z975 + Z80) * r3["se"]     # 80% power 最小可检测效应
    print(f"  [{cut}] n_obs={r3['n_obs']:,} n_ent={r3['n_ent']:,}")
    print(f"     copd main      β={res['copd']['coef']:+.4f} se={res['copd']['se']:.4f} p={res['copd']['p']:.4f}")
    print(f"     copd×age       β={res['c_ag']['coef']:+.4f} p={res['c_ag']['p']:.4f}")
    print(f"     copd×pace      β={res['c_pa']['coef']:+.4f} p={res['c_pa']['p']:.4f}")
    print(f"     copd×age×pace  β={r3['coef']:+.5f} se={r3['se']:.5f} p={r3['p']:.4f}  "
          f"MDβ80={md_beta:+.5f}")
    rows15a.append(dict(cut=cut, term="copd×age×pace", beta=r3["coef"], se=r3["se"],
                        p=r3["p"], ci_lo=r3["ci_lo"], ci_hi=r3["ci_hi"],
                        n_obs=r3["n_obs"], n_ent=r3["n_ent"], MDbeta80=md_beta))
    for t in ["copd","c_ag","c_pa"]:
        rows15a.append(dict(cut=cut, term=t, beta=res[t]["coef"], se=res[t]["se"],
                            p=res[t]["p"], ci_lo=res[t]["ci_lo"], ci_hi=res[t]["ci_hi"],
                            n_obs=res[t]["n_obs"], n_ent=res[t]["n_ent"], MDbeta80=np.nan))

# 4-cell 简单 TWFE（age×pace 二分）作方向性展示
print("  [4-cell simple TWFE for direction]")
cell_rows = []
for gcut, gc in [("age<59","age_grp59"), ("age>=59-inv","age_grp59")]:
    pass
# build 2x2: younger/older x fast/slow pace (median split within pace-merged)
pmed = panel.loc[panel["pace"].notna(), "pace"].median()
panel["pace_fast"] = (panel["pace"] >= pmed).astype(float)
for ag, alab in [(0,"older(>=59)"), (1,"younger(<59)")]:
    for pf, plab in [(0,"slow_pace"), (1,"fast_pace")]:
        dd = panel[(panel["age_grp59"]==ag) & (panel["pace_fast"]==pf)].dropna(subset=[COG,COPD]+COVS).copy()
        dd["copd"]=dd[COPD]
        r = fe_terms(dd, COG, ["copd"], COVS)
        if r:
            print(f"     {alab:12s} {plab:10s} β={r['copd']['coef']:+.4f} p={r['copd']['p']:.4f} n={r['copd']['n_obs']:,}")
            cell_rows.append(dict(stratum=f"{alab}|{plab}", beta=r["copd"]["coef"], p=r["copd"]["p"], n_obs=r["copd"]["n_obs"]))

# ════════════════════════════════════════════════════════════════════════
# Stage 15B : 分年龄层 joint mediation（方向开放）
# ════════════════════════════════════════════════════════════════════════
print("\n=== Stage 15B: age-stratified joint mediation (direction-open) ===")
hcap = make_hhidpn(load_dta(HCAP))   # 4539x13, load full (usecols drops pnfl)
vbs  = make_hhidpn(load_dta(VBS))     # 4611x7,  load full
for src,dst in [("pnfl","nfl"),("pptau181","ptau181")]:
    if src in hcap: hcap[dst]=hcap[src].clip(lower=1e-6); hcap[f"log_{dst}"]=np.log(hcap[dst])
if "phcy" in vbs: vbs["homocys"]=vbs["phcy"]
hcap_c = hcap[["hhidpn","log_nfl","log_ptau181"]].dropna(subset=["log_nfl","log_ptau181"]).drop_duplicates("hhidpn")
vbs_c  = vbs[["hhidpn","homocys"]].dropna(subset=["homocys"]).drop_duplicates("hhidpn")
w13 = panel[panel["wave"]==13].merge(hcap_c, on="hhidpn", how="inner").merge(vbs_c, on="hhidpn", how="inner")
w13 = w13.merge(base.rename("bage"), left_on="hhidpn", right_index=True, how="left")
meds = ["log_nfl","log_ptau181","homocys"]
jmc = w13.dropna(subset=[COG,COPD]+meds+[c for c in COVS if c in w13.columns])
print(f"  joint complete-case n={len(jmc):,}")
covs_m = [c for c in COVS if c in jmc.columns]
rows15b = []
for glab, gmask in [("younger(<59)", jmc["bage"]<59), ("older(>=59)", jmc["bage"]>=59)]:
    dd = jmc[gmask]
    r = joint_mediation(dd, meds, covs_m)
    if r is None:
        print(f"  [{glab}] FAILED"); continue
    print(f"  [{glab}] n={r['n']:,} total={r['total']:+.4f} direct={r['direct']:+.4f}")
    print(f"     JOINT indirect={r['joint']:+.4f} CI[{r['joint_ci'][0]:+.4f},{r['joint_ci'][1]:+.4f}] p={r['joint_p']:.4f}")
    rows15b.append(dict(age_group=glab, pathway="JOINT(all 3)", indirect=r["joint"],
                       ci_lo=r["joint_ci"][0], ci_hi=r["joint_ci"][1], p=r["joint_p"], n=r["n"]))
    for m,lab in [("log_nfl","NfL(log)"),("log_ptau181","pTau181(log)"),("homocys","Homocysteine")]:
        ie,lo,hi,p = r["partial"][m]
        print(f"       {lab:14s} partial={ie:+.4f} CI[{lo:+.4f},{hi:+.4f}] p={p:.4f}")
        rows15b.append(dict(age_group=glab, pathway=lab, indirect=ie, ci_lo=lo, ci_hi=hi, p=p, n=r["n"]))

# ════════════════════════════════════════════════════════════════════════
# Telomere 全样本效应修饰（cog27，无年龄分层）
# ════════════════════════════════════════════════════════════════════════
print("\n=== Telomere full-sample effect modification (cog27, no age stratification) ===")
rows_tel = []
for tlcol, tlname in [("telomere_short", "telomere_short(0/1)"), ("telomere_log", "telomere_log(cont)")]:
    d = panel.dropna(subset=[COG,COPD,tlcol]+COVS).copy()
    d["copd"]=d[COPD]; d["c_tel"]=d["copd"]*d[tlcol]
    res = fe_terms(d, COG, ["copd","c_tel"], COVS)
    if res is None:
        print(f"  [{tlname}] FAILED"); continue
    rt = res["c_tel"]
    print(f"  [{tlname}] copd×tel β={rt['coef']:+.4f} se={rt['se']:.4f} p={rt['p']:.4f} n={rt['n_obs']:,}")
    rows_tel.append(dict(modifier=tlname, term="copd×telomere", beta=rt["coef"], se=rt["se"],
                         p=rt["p"], ci_lo=rt["ci_lo"], ci_hi=rt["ci_hi"], n_obs=rt["n_obs"]))

# ════════════════════════════════════════════════════════════════════════
print("\nSAVE")
pd.DataFrame(rows15a).to_csv(OUT/"stage15A_threeway.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(rows15b).to_csv(OUT/"stage15B_age_mediation.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(rows_tel).to_csv(OUT/"telomere_full_modification.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(cell_rows).to_csv(OUT/"stage15A_4cell.csv", index=False, encoding="utf-8-sig")
print("DONE")
