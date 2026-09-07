# -*- coding: utf-8 -*-
"""
STEP 1 (优先级最高) — 决定性 cog27 口径 59 岁分半 TWFE + 组间 Wald
====================================================================
决定标题强度: "reversal"(显著) vs "age-dependent heterogeneity"(边缘/不显著)。
附带 09C vs Stage14 older 组估计矛盾的规格敏感性网格。

数据: Stage2 面板 wave 上限=13(无 14/15)，"全样本"=waves 3-13。
      主结局 cog27(148,760 非空)；Stage14 用 cogtot(114,222)。
      baseline age = 每人 waves 3-13 首次观测 agey_b。
      主切分点 = 59(稿件 Table2 脚注 + 用户指定)；实测中位数=61，作为敏感性另报。
模型: 双向 demean(hhidpn+wave) + hhidpn 聚类稳健 SE(同 Stage12/13/14)。
"""
import numpy as np, pandas as pd, time
from scipy import stats

PANEL = r"D:\HRSCOPD\数据分析\COPD与认知\投稿代码\code\Pipeline_Config_Refactor\输出结果\hrs_panel_stage2.parquet"
OUT_CSV = r"D:\HRSCOPD\数据分析\COPD与认知\方案A\step1_cog27_reversal.csv"

def twoway_demean(df, y_col, x_cols, entity="hhidpn", time="wave", max_iter=300, tol=1e-9):
    cols = [entity, time, y_col] + x_cols
    d = df[[c for c in cols if c in df.columns]].copy()
    for c in [y_col] + x_cols:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna()
    if len(d) < 50: return None
    eid, tid = d[entity].values, d[time].values
    Y = d[y_col].values.astype(float); X = d[x_cols].values.astype(float)
    ec, el = pd.factorize(eid); tc, tl = pd.factorize(tid)
    ne, nt = len(el), len(tl)
    def _dm(v):
        v = v - v.mean()
        for _ in range(max_iter):
            v0 = v.copy()
            v = v - (np.bincount(ec, weights=v, minlength=ne)/np.maximum(np.bincount(ec, minlength=ne),1))[ec]
            v = v - (np.bincount(tc, weights=v, minlength=nt)/np.maximum(np.bincount(tc, minlength=nt),1))[tc]
            if np.max(np.abs(v-v0)) < tol: break
        return v
    return _dm(Y.copy()), np.column_stack([_dm(X[:,j].copy()) for j in range(X.shape[1])]), eid, len(d)

def clustered_se(Xd, resid, clusters):
    n, k = Xd.shape; XtXi = np.linalg.pinv(Xd.T@Xd); G = len(np.unique(clusters))
    meat = np.zeros((k,k))
    for cl in np.unique(clusters):
        m = clusters==cl; s = Xd[m].T@resid[m]; meat += np.outer(s,s)
    corr = (G/(G-1))*((n-1)/(n-k))
    diag = np.diag(corr*XtXi@meat@XtXi)
    return np.sqrt(np.maximum(diag,0))   # 防御数值负对角

def run_fe(df, outcome, x_cols, entity="hhidpn", time="wave"):
    """返回 {var:(coef,se,p)} 与 n_obs/n_ent"""
    res = twoway_demean(df, outcome, x_cols, entity, time)
    if res is None: return None
    Yd, Xd, eid, n = res; k = Xd.shape[1]
    if n < 100 or n <= k: return None
    betas = np.linalg.lstsq(Xd, Yd, rcond=None)[0]
    resid = Yd - Xd@betas; n_ent = len(np.unique(eid))
    ses = clustered_se(Xd, resid, eid)
    out = {"n_obs": n, "n_ent": n_ent}
    for i, name in enumerate(x_cols):
        coef, se = betas[i], ses[i]
        t = coef/se if se>0 and not np.isnan(se) else np.nan
        p = 2*(1-stats.t.cdf(abs(t), max(n-n_ent-k,1))) if not np.isnan(t) else np.nan
        out[name] = (coef, se, p)
    return out

def stars(p): return "***" if (p is not None and p<.001) else "**" if (p is not None and p<.01) else "*" if (p is not None and p<.05) else "n.s."

def wald(c1,s1,c2,s2):
    d=c1-c2; se=np.sqrt(s1**2+s2**2); z=d/se if se>0 else np.nan
    p=2*(1-stats.norm.cdf(abs(z))) if not np.isnan(z) else np.nan
    return d,z,p

t0=time.time()
panel = pd.read_parquet(PANEL)
print(f"panel {panel.shape} persons {panel['hhidpn'].nunique():,} waves {panel['wave'].min()}-{panel['wave'].max()}")
print(f"  cog27 nonnull={panel['cog27'].notna().sum():,}  cogtot nonnull={panel['cogtot'].notna().sum():,}")

FULL=["agey_b","female","raedyrs","race_white","race_black","hispanic","married","smoken","cesd","bmi","hibpe","diabe","stroke"]
MIN=["bmi","conde","hibpe","diabe","hearte","stroke","cancre","smoken","shlt"]
FULL=[c for c in FULL if c in panel.columns]; MIN=[c for c in MIN if c in panel.columns]

sub=panel[panel["wave"].between(3,13)].sort_values(["hhidpn","wave"])
bl=sub.groupby("hhidpn")["agey_b"].first().rename("baseline_age")
panel=panel.merge(bl,on="hhidpn",how="left")
emp_med=panel["baseline_age"].median()
print(f"  实测 baseline-age 中位数 = {emp_med:.2f}  (稿件/用户指定切分点 = 59)")

panel["older59"]=(panel["baseline_age"]>=59).astype(float)
panel["older_emp"]=(panel["baseline_age"]>=emp_med).astype(float)
panel["copd_x_59"]=panel["copd_ever"]*panel["older59"]
panel["copd_x_emp"]=panel["copd_ever"]*panel["older_emp"]
panel["age_c59"]=(panel["baseline_age"]-59).astype(float)
panel["age_cemp"]=(panel["baseline_age"]-emp_med).astype(float)
panel["copd_x_age59"]=panel["copd_ever"]*panel["age_c59"]
panel["copd_x_ageemp"]=panel["copd_ever"]*panel["age_cemp"]

rows=[]
def report_split(cog, spec, covs, split_flag, xinter):
    psub=panel[panel["wave"].between(3,13)]
    ry=run_fe(psub[psub[split_flag]==0], cog, ["copd_ever"]+covs)
    ro=run_fe(psub[psub[split_flag]==1], cog, ["copd_ever"]+covs)
    if not ry or not ro:
        print(f"  [{cog}|{spec}] 分半样本不足，跳过"); return
    cy,sy,py=ry["copd_ever"]; co,so,po=ro["copd_ever"]
    d,z,p=wald(cy,sy,co,so)
    rows.append(dict(outcome=cog,spec=spec,split="59" if "59" in split_flag else "emp",
                     model="split-Wald",younger_beta=cy,younger_se=sy,younger_p=py,
                     older_beta=co,older_se=so,older_p=po,wald_delta=d,wald_z=z,wald_p=p,
                     n_young=ry["n_obs"],n_older=ro["n_obs"]))
    print(f"\n[{cog} | {spec} | split@{('59' if '59' in split_flag else 'emp')}]")
    print(f"  younger β={cy:+.4f} SE={sy:.4f} p={py:.4f} {stars(py)} n={ry['n_obs']:,}")
    print(f"  older   β={co:+.4f} SE={so:.4f} p={po:.4f} {stars(po)} n={ro['n_obs']:,}")
    print(f"  WALD Δ={d:+.4f} z={z:+.3f} p={p:.4f} -> {'REVERSAL(显著)' if (p is not None and p<.05) else 'age-heterogeneity(边缘/不显著)'}")
    rp=run_fe(psub, cog, ["copd_ever",split_flag,xinter]+covs)
    if rp:
        cb,sb,pb=rp[xinter]
        rows.append(dict(outcome=cog,spec=spec,split="59" if "59" in split_flag else "emp",
                         model="pooled-interaction",interaction_beta=cb,interaction_se=sb,interaction_p=pb,n_obs=rp["n_obs"]))
        print(f"  pooled 交互(coef on {xinter}) β={cb:+.4f} SE={sb:.4f} p={pb:.4f} {stars(pb)}")

print("\n############ 主分析: cog27, 切分点=59 (稿件口径) ############")
for spec,covs in [("FULL(incl agey_b)",FULL),("MIN(09C,no agey_b)",MIN)]:
    report_split("cog27", spec, covs, "older59", "copd_x_59")
print("\n############ cog27 敏感性: 实测中位数切分 ############")
for spec,covs in [("FULL(incl agey_b)",FULL),("MIN(09C,no agey_b)",MIN)]:
    report_split("cog27", spec, covs, "older_emp", "copd_x_emp")
print("\n############ 交叉验证: cogtot (Stage14 同款结局) ############")
for spec,covs in [("FULL(incl agey_b)",FULL),("MIN(09C,no agey_b)",MIN)]:
    report_split("cogtot", spec, covs, "older59", "copd_x_59")
print("\n############ 连续年龄梯度(稳健性): cog27 FULL ############")
psub=panel[panel["wave"].between(3,13)]
for lbl,xint in [("59-centered","copd_x_age59"),("emp-centered","copd_x_ageemp")]:
    rc=run_fe(psub,"cog27",["copd_ever","age_c59" if "59" in xint else "age_cemp",xint]+FULL)
    if rc:
        cb,sb,pb=rc[xint]
        rows.append(dict(outcome="cog27",spec="FULL",split=lbl,model="continuous-age-grad",
                         interaction_beta=cb,interaction_se=sb,interaction_p=pb,n_obs=rc["n_obs"]))
        print(f"  [cog27 FULL {lbl}] copd×(age-59c) β={cb:+.4f} SE={sb:.4f} p={pb:.4f} {stars(pb)}")

df=pd.DataFrame(rows); df.to_csv(OUT_CSV,index=False,encoding="utf-8-sig")
print(f"\nSaved {OUT_CSV} ({len(df)} rows) in {time.time()-t0:.1f}s")
