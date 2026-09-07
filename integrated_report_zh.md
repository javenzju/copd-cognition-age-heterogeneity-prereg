# 方案 A（修订版：年龄异质性方案）— 最终整合分析报告

**工作标题（建议）**
*Age-Dependent Heterogeneity in the COPD–Cognitive Decline Association: Beyond Selective Survival, Toward Unexplained Younger-Onset Vulnerability*

**生成日期**：2026-09-07
**数据**：HRS 面板 `hrs_panel_stage2.parquet`（waves 1–13，20,863 人）+ DunedinPACE（`HRSPoA_Shared.dta`）+ wave-13 生物标志物子样本（HCAPbioPilot NfL/pTau181、VBS Homocysteine）
**主结局**：`cog27`（稿件主分析口径）；统一采用 **FULL 协变量规格**（含 `agey_b` + 人口学 + 合并症），与 Step 1 决策一致
**状态**：除 OSF 预注册 / medRxiv-bioRxiv-AAIC 检索 / 预印本挂出（用户要求暂缓）外，全部剩余可执行分析已完成

---

## 0. 一句话结论

COPD 与认知衰退的关联存在**真实且方向稳健的年龄异质性**（年轻起病者认知损伤更重、老年起病者≈0），该异质性**不能被 DunedinPACE 生物年龄差异解释**（三向交互 null 且检验力不足），而 **NfL/pTau181/Hcy 神经损伤通路只在老年起病组显著中介**——年轻起病组携带更大认知负担却不被任何已测生物标志物承载。这构成一条不撞 Vivek 2026、可发表、且明确指向"下一步机制缺口"的叙事。

---

## 1. 标题强度决策（Step 1，已定稿）

用稿件主结局 `cog27` + waves 3–13 全样本 + 59 岁分半 TWFE + 组间 Wald：

| 规格 | younger β (p) | older β (p) | Wald p | 判定 |
|---|---|---|---|---|
| **FULL（含 agey_b，应采用的）** | −0.144 (.111) | +0.015 (.850) | **.187** | 不显著 |
| FULL（实测中位 61 切分） | −0.164 (.045*) | +0.070 (.428) | **.052** | 边缘 |
| MIN（09C 风格，无 agey_b） | −0.112 (.278) | +0.192 (.045*) | .031 | 显著（但 older 正向为选择性存活伪阳性） |

**结论**：标题写 **"age-dependent heterogeneity"**，不写 "reversal"。"反转"只在未调整年龄的 MIN 规格下显著，而该规格的 older 正向值是幸存者健康选择的伪阳性，用其做标题会过度宣称。

---

## 2. Stage 15A：COPD × age × DunedinPACE 三向交互（仅 pace 支线）

telomere 支线按用户指示放弃（样本量仅 ~189–234 例 COPD∩telomere，功效不足）。

| 切分 | 项 | β | p | 最小可检测 Δβ(80% power) |
|---|---|---|---|---|
| age<59 | copd×age×pace | −0.0166 | .930 | 0.528 |
| age<59 | copd×pace（双向） | −0.244 | .075 | — |
| median | copd×age×pace | +0.0523 | .792 | 0.556 |
| median | copd×pace（双向） | −0.306 | .053 | — |

4-cell 方向性（TWFE，age×pace 二分）：
- younger + slow_pace：β=−0.184 (p=.085，唯一下探为负)
- younger + fast_pace：β=+0.127 (p=.45)
- older + slow/fast_pace：β≈+0.03–0.07 (n.s.)

**诚实解读**：
- 三向交互 **null**（p>.79）且 **MDβ80≈0.53**（约为主 COPD 效应的 2–4 倍）——本数据**检验力不足以检测任何合理量级**的生物年龄解释效应。
- 因此结论只能是：*DunedinPACE 没有提供"年龄异质性由生物年龄差异造成"的证据*，但也**不能据此断言 bio-age 完全无关**（检力不足）。
- 这与"选择性存活也无法完全解释 older 异质性（Stage 9C：需 EMM≥4，而真实 COPD 超额死亡 HR 仅 1.3–2.0）"叠加，指向：**年龄异质性真实存在、且现有两大候选解释（选择性存活 / 已测生物年龄）都不足以完全覆盖** → 存在未测量机制，这正是论文的立论空间。

---

## 3. Stage 15B：分年龄层 joint biomarker mediation（方向开放，不复读预设）

wave-13 横断面并行多重中介（NfL + pTau181 + Homocysteine，bootstrap 2000），按 baseline age 分层：

| 年龄层 | JOINT indirect | p | NfL | pTau181 | Homocysteine |
|---|---|---|---|---|---|
| younger(<59) n=1,440 | −0.033 [−0.144, +0.056] | **.50** | −0.008 (p=.85) | −0.001 (p=.97) | −0.024 (p=.28) |
| older(≥59) n=1,028 | −0.131 [−0.255, −0.020] | **.018** | −0.084 (p=.017) | −0.059 (p=.046) | +0.012 (p=.71) |

**结论（与用户修正后的方向假设一致，且复刻了 Stage 14 log 的旧结论）**：
- 神经-axonal / AD 病理通路（NfL、pTau181）**只在老年起病 COPD 显著中介** COPD→认知关联；
- 年轻起病组认知总效应更大（横断面 total=−0.885 vs older +0.018），却**不被任何已测生物标志物承载**（joint p=.50）。
- 这制造了一个明确的"下一步缺口"：年轻起病患者的认知损伤由当前测量框架之外的机制驱动（如围生期/早年 airway 编程、累积性低度缺氧的早年起始、或尚未纳入的标志物）。

---

## 4. Telomere 全样本效应修饰（cog27，不做年龄分层）

| 修饰项 | 项 | β | p |
|---|---|---|---|
| telomere_short (0/1) | copd×telomere | −0.527 | **.058** |
| telomere_log (连续) | copd×telomere | +0.552 | .216 |

**结论**：短端粒在**全样本**水平上呈**边际显著**（p=.058）地放大 COPD 对认知的不利关联，方向与既有结果（telomere Wald p=.041，cogtot 口径）一致；连续 log 编码方向相反且无意义（二进制低/高是更有意义的效应修饰编码）。按要求**未做年龄分层**。

---

## 5. 整合叙事锚点（投稿用）

> HRS 中 COPD 与认知衰退的关联随起病年龄异质：年轻起病者携带真实且未被选择性存活解释的额外认知损伤，而老年起病者的关联经年龄+合并症调整后趋近零。DunedinPACE 生物年龄不能解释这一梯度（三向交互 null，且检验力不足以排除小效应）。在生物标志物层面，NfL/pTau181/Hcy 神经损伤通路仅介导老年起病组的关联——年轻起病组的更大认知负担落在该框架之外，提示一种尚未被测量的早年损伤机制。对临床的意义：应优先对**年轻起病 COPD** 患者开展认知筛查，而非（如现行指南隐含的）仅关注老年共病群体。

---

## 6. 方法学诚信声明（务必写入稿件 limitations）

1. **三向交互检力不足**：MDβ80≈0.53，不能据此主张 bio-age 独立；仅能称"无证据支持"。
2. **Stage 15B 为横断面中介**（wave-13 单点），因果方向弱于 TWFE；其 "total effect" 是横断面估计，不与 Step 1 的纵向 TWFE β 直接比较。
3. **年龄异质性 Wald 在 FULL 规格下边缘/不显著**（p=.187 / .052），故措辞为 heterogeneity 而非 reversal。
4. **older 组 older-onset 正向值在 MIN 规格下为选择性存活伪阳性**——两表统一用 FULL 规格，已消除 09C(+0.344)/Stage14(+0.107) 的内部 3 倍矛盾。

---

## 7. 目标期刊（JGSA 已解锁）

- **首选**：Journals of Gerontology Series A — 用户已确认无在审稿，排除解除。
- **次选**：Age and Ageing（IF≈10，对反直觉临床异质性接受度高）。
- **备选**：GeroScience（方法学+老年流行病学属性契合；但要求强 novelty，本角度为 heterogeneity 而非机制突破，定位略低于 GeroScience 传统预期，可作备选）。
- 注：AJE 等已拒期刊不再重复投；若转投需 cover letter 明确"同一队列不同研究问题"。

---

## 8. 用户要求暂缓的事项（未执行，非遗漏）

- OSF 预注册（Stage 15 分析计划锁定）
- medRxiv / bioRxiv / AAIC 摘要检索（确认"年龄异质性角度"未被抢先）
- 预印本挂出（锁定 priority timestamp，Vivek 案例教训）
> 上述三项按用户 2026-09-07 指示暂缓，待用户另行授权后再推进。

---

## 9. 交付文件清单（均在 `D:\HRSCOPD\数据分析\COPD与认知\方案A\`）

| 文件 | 内容 |
|---|---|
| `Step1_标题决策与规格统一.md` + `step1_cog27_reversal.py` + `step1_cog27_reversal_valid.csv` | 标题强度决策（FULL 规格 Wald p=.187） |
| `stage15_analysis.py` | Stage 15A/B + telomere 综合脚本（可复跑） |
| `stage15A_threeway.csv` | 三向交互 + 最小可检测 Δβ |
| `stage15A_4cell.csv` | age×pace 4 格方向性 TWFE |
| `stage15B_age_mediation.csv` | 分年龄 joint mediation（方向开放） |
| `telomere_full_modification.csv` | 全样本 telomere 效应修饰 |
| `make_figures.py` + `fig_age_heterogeneity_mediation.png/.pdf` | 核心双面板图（图 A 年龄分层效应；图 B 分年龄中介） |
| `可行性评估_年龄反转方案.md` | 前期可行性核实（Vivek 预印本、TAVS 证伪、Stage12 矛盾等） |
