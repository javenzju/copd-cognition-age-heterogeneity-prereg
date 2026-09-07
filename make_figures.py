# -*- coding: utf-8 -*-
"""方案A 核心两张图：年龄分层 COPD 效应 + 分年龄 joint mediation。
读取 step1 / stage15B 的 CSV，生成 publication-ready 双面板图。"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
OUT = Path(r"D:\HRSCOPD\数据分析\COPD与认知\方案A")

def ci_from_bp(beta, p, df=1e6):
    if p <= 0 or p >= 1: return (np.nan, np.nan)
    t = stats.t.ppf(1 - p/2, df)
    se = abs(beta)/t
    return (beta-1.96*se, beta+1.96*se)

s1 = pd.read_csv(OUT/"step1_cog27_reversal_valid.csv")
row = s1[(s1.outcome=="cog27") & (s1.spec.str.startswith("FULL")) & (s1.split=="59")].iloc[0]
y_b, y_p = row.younger_beta, row.younger_p
o_b, o_p = row.older_beta,   row.older_p
y_ci = ci_from_bp(y_b, y_p); o_ci = ci_from_bp(o_b, o_p)

m = pd.read_csv(OUT/"stage15B_age_mediation.csv")
yj = m[(m.age_group=="younger(<59)") & (m.pathway=="JOINT(all 3)")].iloc[0]
oj = m[(m.age_group=="older(>=59)")  & (m.pathway=="JOINT(all 3)")].iloc[0]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.6))
fig.patch.set_facecolor("white")

# Panel A: age-stratified TWFE effect
axA.axvline(0, color="#888", ls="--", lw=1)
for i,(lab,b,ci,p) in enumerate([("Younger-onset (baseline age <59)", y_b, y_ci, y_p),
                                 ("Older-onset (baseline age ≥59)",  o_b, o_ci, o_p)]):
    col = "#1b4f72"
    axA.plot([ci[0],ci[1]], [i,i], color=col, lw=2, solid_capstyle="round")
    axA.plot(b, i, "o", color=col, ms=9,
             markerfacecolor=("white" if p>=0.05 else col), markeredgewidth=1.8)
    sig = "n.s." if p>=0.05 else ("*" if p<0.05 else "**")
    axA.text(ci[1]+0.02, i, f"β={b:+.3f} [{ci[0]:+.3f}, {ci[1]:+.3f}] p={p:.3f} {sig}",
             va="center", fontsize=8.5, color=col)
axA.set_yticks([0,1]); axA.set_yticklabels(
    ["Younger-onset\n(baseline age <59)","Older-onset\n(baseline age ≥59)"], fontsize=9)
axA.set_xlabel("TWFE β (COPD × time → cog27), 95% CI", fontsize=9)
axA.set_xlim(-0.45, 0.55)
axA.set_title("Panel A. Age-dependent heterogeneity in the COPD–cognition association",
              fontsize=9.5, fontweight="bold", loc="left")
axA.spines[["top","right"]].set_visible(False)
axA.text(0.0, -0.18, f"Between-group Wald Δβ={row.wald_delta:+.3f}, p={row.wald_p:.3f}  "
         f"(FULL specification, outcome = cog27)",
         transform=axA.transAxes, fontsize=7.5, color="#555")

# Panel B: joint mediation indirect effect by age
axB.axvline(0, color="#888", ls="--", lw=1)
rowsB = [("Younger-onset", yj.indirect, yj.ci_lo, yj.ci_hi, yj.p),
         ("Older-onset",   oj.indirect, oj.ci_lo, oj.ci_hi, oj.p)]
for i,(lab,b,lo,hi,p) in enumerate(rowsB):
    col = "#b03a2e" if (i==1 and p<0.05) else "#1b4f72"
    axB.plot([lo,hi], [i,i], color=col, lw=2, solid_capstyle="round")
    axB.plot(b, i, "o", color=col, ms=9,
             markerfacecolor=("white" if p>=0.05 else col), markeredgewidth=1.8)
    sig = "n.s." if p>=0.05 else ("*" if p<0.05 else "**")
    axB.text(hi+0.01, i, f"ab={b:+.3f} [{lo:+.3f}, {hi:+.3f}] p={p:.3f} {sig}",
             va="center", fontsize=8.5, color=col)
axB.set_yticks([0,1]); axB.set_yticklabels(["Younger-onset\n(baseline age <59)","Older-onset\n(baseline age ≥59)"], fontsize=9)
axB.set_xlabel("Joint indirect effect (NfL + pTau181 + Homocysteine), 95% bootstrap CI", fontsize=9)
axB.set_xlim(-0.30, 0.20)
axB.set_title("Panel B. Neuro-biomarker mediation is significant only in older-onset COPD",
              fontsize=9.5, fontweight="bold", loc="left")
axB.spines[["top","right"]].set_visible(False)
axB.text(0.0, -0.18,
         f"Joint model n={int(yj.n)} (younger) / {int(oj.n)} (older); bootstrap 2000; "
         f"younger JOINT p={yj.p:.2f} vs older JOINT p={oj.p:.3f}",
         transform=axB.transAxes, fontsize=7.5, color="#555")

plt.tight_layout(rect=[0,0.04,1,1])
fig.savefig(OUT/"fig_age_heterogeneity_mediation.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT/"fig_age_heterogeneity_mediation.pdf", bbox_inches="tight")
print("fig saved")
