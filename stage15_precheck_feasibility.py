# -*- coding: utf-8 -*-
"""
Stage-15 预检可行性表 (等价物)
================================
回答: 1) telomere / DunedinPACE 在两年龄层的分布是否不同(生物年龄解释有无戏)
      2) 两年龄层内 COPD 病例数、telomere/pace 重叠 cell 数 → 三向交互功效底线
      3) telomere-pace 相关(复验 Stage12 r=-0.024)
口径: baseline age = wave3-15 首条 agey_b; split at 59 (Stage14/09C 共用) 和 中位数 两版
      COPD = copd_ever (ever); 附 2016 横断面 COPD (wave13) 计数
"""
import numpy as np, pandas as pd
from scipy import stats

PANEL = r"D:\HRSCOPD\数据分析\COPD与认知\投稿代码\code\Pipeline_Config_Refactor\输出结果\hrs_panel_stage2.parquet"
POA   = r"D:\HRS\PaceOfAging\HRSPoA_Shared.dta"

panel = pd.read_parquet(PANEL)
poa, _ = __import__("pyreadstat").read_dta(POA)
poa["hhidpn"] = pd.to_numeric(poa["hhidpn"], errors="coerce").astype("Int64")
poa = poa[["hhidpn", "paceofaging"]].dropna()
poa["paceofaging"] = poa["paceofaging"].clip(poa["paceofaging"].quantile(0.01), poa["paceofaging"].quantile(0.99))

# baseline age (wave3-15 first)
sub = panel[panel["wave"].between(3, 15)].sort_values(["hhidpn", "wave"])
bl = sub.groupby("hhidpn")["agey_b"].first().rename("baseline_age")
panel = panel.merge(bl, on="hhidpn", how="left")

# person-level key vars from panel
p = panel.drop_duplicates("hhidpn")
p["tel_ts"] = p.groupby("hhidpn")["telomere_ts"].transform("first")
p["copd_ever"] = p.groupby("hhidpn")["copd_ever"].transform("max")
p["copd_2016"] = p.groupby("hhidpn")["copd_wave"].transform(lambda s: 1 if (s == 13).any() else 0)
# 更接近稿件定义: wave13 copd? 用 copd_wave==13
pp = p[["hhidpn", "baseline_age", "tel_ts", "copd_ever"]].drop_duplicates("hhidpn")
pp = pp.merge(poa, on="hhidpn", how="left")
pp["w13"] = panel[panel["wave"] == 13].drop_duplicates("hhidpn").set_index("hhidpn")["copd_wave"].reindex(pp["hhidpn"]).values
pp["copd_w13"] = (pp["w13"] == 13).astype(int)

def cross(pp, splitval, label):
    pp = pp.copy()
    pp["grp"] = np.where(pp["baseline_age"] < splitval, "younger(<%d)" % splitval, "older(>=%d)" % splitval)
    print(f"\n########## split = {label}: baseline age < {splitval} ##########")
    print(f"  baseline-age median: {pp['baseline_age'].median():.1f}")
    for g, d in pp.groupby("grp"):
        n = len(d)
        copd = int(d["copd_ever"].sum())
        copd13 = int(d["copd_w13"].sum())
        tel = int(d["tel_ts"].notna().sum())
        pace = int(d["paceofaging"].notna().sum())
        telcopd = int((d["tel_ts"].notna() & (d["copd_ever"] == 1)).sum())
        pacecopd = int((d["paceofaging"].notna() & (d["copd_ever"] == 1)).sum())
        both = int((d["tel_ts"].notna() & d["paceofaging"].notna()).sum())
        print(f"  {g:14s} n={n:6,}  COPD_ever={copd:5,}  COPD(w13)={copd13:4,}  "
              f"telomere={tel:6,}  pace={pace:6,}  tel&copd={telcopd:4,}  pace&copd={pacecopd:4,}  tel&pace={both:5,}")
        sub = d[d["grp"] == g]
        if sub["tel_ts"].notna().sum() > 30:
            print(f"      telomere T/S mean={sub['tel_ts'].mean():.3f} SD={sub['tel_ts'].std():.3f}")
        if sub["paceofaging"].notna().sum() > 30:
            print(f"      DunedinPACE mean={sub['paceofaging'].mean():.3f} SD={sub['paceofaging'].std():.3f}")
    return pp

pp = cross(pp, 59, "Stage14/09C 口径")
# 再按实际中位数
med = int(pp["baseline_age"].median())
pp = cross(pp, med, "数据中位数")

# telomere-pace correlation overall & by split
d = pp[["tel_ts", "paceofaging"]].dropna()
if len(d) > 100:
    r, pv = stats.pearsonr(d["tel_ts"], d["paceofaging"])
    print(f"\ntelomere x DunedinPACE: r={r:+.4f}  p={pv:.4f}  n={len(d):,}  (Stage12报 r=-0.024 p=.207)")
d = pp[["baseline_age", "tel_ts"]].dropna(); r, pv = stats.pearsonr(d["baseline_age"], d["tel_ts"])
print(f"telomere x baseline_age: r={r:+.4f}  p={pv:.4f}  n={len(d):,}")
d = pp[["baseline_age", "paceofaging"]].dropna(); r, pv = stats.pearsonr(d["baseline_age"], d["paceofaging"])
print(f"pace x baseline_age:     r={r:+.4f}  p={pv:.4f}  n={len(d):,}")

# telomere / pace by COPD within split (59)
print("\n########## 按 COPD × age 的生物年龄分布 (split=59) ##########")
pp["grp"] = np.where(pp["baseline_age"] < 59, "younger", "older")
for g in ["younger", "older"]:
    for c in [0, 1]:
        sub = pp[(pp["grp"] == g) & (pp["copd_ever"] == c)]
        line = f"  {g:8s} COPD={c}: n={len(sub):6,}"
        if sub["tel_ts"].notna().sum() > 30:
            line += f"  tel mean={sub['tel_ts'].mean():.3f}"
        if sub["paceofaging"].notna().sum() > 30:
            line += f"  pace mean={sub['paceofaging'].mean():.3f}"
        print(line)
# COPD status vs bioage within age group (t-test, adjusted later; simple here)
for g in ["younger", "older"]:
    for var in ["tel_ts", "paceofaging"]:
        a = pp[(pp["grp"] == g) & (pp["copd_ever"] == 1)][var].dropna()
        b = pp[(pp["grp"] == g) & (pp["copd_ever"] == 0)][var].dropna()
        if len(a) > 20 and len(b) > 20:
            t, pv = stats.ttest_ind(a, b, equal_var=False)
            print(f"  {g:8s} {var:12s}: COPD={len(a):4d}(mean {a.mean():.3f}) vs non={len(b):6d}(mean {b.mean():.3f})  t={t:+.2f} p={pv:.4f}")
