import pandas as pd, numpy as np
from scipy.stats import beta, spearmanr, wilcoxon
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from pathlib import Path
ROOT = Path("/Users/donny/LOY_UM")  # local clone of this repository; edit if cloned elsewhere
DATA, RESULTS, FIGS = ROOT / "data", ROOT / "results", ROOT / "figures"
for _d in (RESULTS / "Source_Data", FIGS / "main", FIGS / "supplementary"): _d.mkdir(parents=True, exist_ok=True)
F = RESULTS / "Source_Data" / "LOY_UM_MasterData.xlsx"
x = pd.read_excel(F, None)
s1 = pd.read_csv(DATA / "FigS1_summary.csv")
seg = pd.read_csv(DATA / "seg_027_051.csv").rename(columns={"LX27_pre": "USC-027_tumor", "LX51_pre": "USC-051_tumor"})

# --- README updates
r = x["README"]
r.loc[r.Item == "sWGS chrY window", "Definition"] = "Male-specific region, chrY 2.8–26.6 Mb (coordinates as in the segment files; confirm reference build — Methods state hg38)"
r.loc[r.Item == "Arm calls", "Definition"] = "TFx-adjusted arm mean < 0.75 = loss, > 1.25 = gain (centromere coordinates matching the segment files)"
r.loc[r.Item == "Excluded tumors", "Definition"] = ("USC-027, USC-051: ploidy-indeterminate hyperdiploid profiles (gains only, no losses); tumor chrY:chrX 1.74 and 1.10 "
                                                    "(no relative Y loss); blood retained; 027 blood = positive control (73.9% sWGS, 71.6% dPCR; age 97)")
r = pd.concat([r, pd.DataFrame([
    ["Manuscript version", "LOY_UM_Manuscript_Cancers_Xu_v5 (Figures 1–4; Figures S1–S4; Tables S1–S3; Data S1)"],
    ["Wilcoxon test", "Two-sided exact signed-rank on paired tumor vs blood raw LOY fractions in the 7 LOY-positive patients"],
    ["dPCR threshold meaning", "Analytical calling threshold for the assay, not a biological cutoff for mLOY"]], columns=r.columns)], ignore_index=True)

# --- Index
idx = pd.DataFrame([
    ["Figure 1A", "Fig1_cohort", "Cohort derivation counts"],
    ["Figure 2A", "Fig2A_frequency", "Frequency of LOY and arm-level SCNAs (n = 9)"],
    ["Figure 2B", "Fig2B_genomewide", "Genome-wide TFx-adjusted tumor-cell copy state"],
    ["Figure 2C", "Fig2C_LOY_TFx", "Raw and TFx-adjusted LOY fractions, TFx scenarios"],
    ["Figure 2 / Table S2", "Fig2_arm_copy_state", "Arm-level TFx-adjusted copy states"],
    ["Figure 3A", "Fig3A_tumor_vs_blood", "Paired tumor vs matched blood raw LOY (9 pairs)"],
    ["Figure 3B", "Fig3B_calibration; Fig3B_parameters", "11 matched bloods, baseline, tolerance limits, 027 positive control"],
    ["Figure 3C", "Fig3C_threshold_sweep", "Calls across thresholds 0–30%"],
    ["Figure 4A–B", "Fig4AB_dPCR_cohort; Fig4B_stats", "dPCR mLOY in 61 men, CIs, age"],
    ["Figure 4C", "Fig4C_concordance", "sWGS vs dPCR in 11 bloods"],
    ["Figure 4D", "Fig4D_groups", "Tumor, matched blood, cohort distributions"],
    ["Figure S1", "FigS1_summary; FigS1_seg_027_051", "Excluded ploidy-indeterminate tumors: X/Y dosage and segment data"],
    ["Figure S2A", "FigS2A_blood_baselines", "Matched blood LOY under three baselines"],
    ["Figure S2B", "FigS2B_blood_genomewide", "Genome-wide matched blood profiles"],
    ["Figure S3", "FigS3_agreement; FigS3_parameters", "Bland–Altman and Y:X vs Y:autosome"],
    ["Figure S4", "FigS4_historical", "Historical karyotype re-plot (Prescher 1990, Horsman 1993)"],
    ["Text", "Key_numbers", "Every number quoted in the Abstract/Results, recomputed from source"],
    ["Controls", "dPCR_female_controls", "Female assay-specificity controls"],
    ["Summaries", "blood_sWGS_summary; tumor_sWGS_summary", "Per-sample chrY, chrX, autosome means and LOY fractions"]],
    columns=["Figure/panel", "Sheet(s)", "Content"])

# --- Fig1 cohort
fig1 = pd.DataFrame([
    ["Intraocular specimens screened", 35, ""], ["Female (excluded)", 18, ""], ["Male specimens", 17, "6 TW, 11 AH"],
    ["Indeterminate (no detectable SCNA)", 6, "2 TW, 4 AH"], ["SCNA-positive", 11, ""],
    ["Ploidy-indeterminate (excluded from tumor analysis)", 2, "USC-027, USC-051"], ["Evaluable tumors", 9, "all with matched blood sWGS"],
    ["Expanded blood cohort (dPCR)", 61, "includes all 11 SCNA-positive patients"]], columns=["Step", "n", "Note"])

# --- Key numbers recomputed
t = x["Fig2C_LOY_TFx"]; pr = x["Fig3A_tumor_vs_blood"]; co = x["Fig4AB_dPCR_cohort"]; cc = x["Fig4C_concordance"]; sa = x["FigS3_agreement"]
pos = t[t.LOY_call == "LOY+"]
ci = lambda k, n: (beta.ppf(.025, k, n - k + 1) * 100, beta.ppf(.975, k + 1, n - k) * 100)
k1 = ci(7, 9); k2 = ci(4, 61)
pp = pr[pr.tumor_call == "LOY+"]
w = wilcoxon(pp.tumor_LOY_pct, pp.matched_blood_LOY_pct, method="exact")
rho = spearmanr(co.age_at_diagnosis, co.mLOY_pct); no27 = co[co.Patient != "USC-027"]; rho2 = spearmanr(no27.age_at_diagnosis, no27.mLOY_pct)
ev = cc[~cc.Patient.isin(["USC-027", "USC-051"])]
hw = ((co.mLOY_CI95_high_pct - co.mLOY_CI95_low_pct) / 2).median()
K = [
    ["Tumor LOY-positive", "7/9", f"77.8% (95% CI {k1[0]:.1f}–{k1[1]:.1f}%)", "Fig2C_LOY_TFx"],
    ["Raw LOY, LOY-positive tumors", f"{pos.raw_LOY_fraction_pct.min():.1f}–{pos.raw_LOY_fraction_pct.max():.1f}%", "", "Fig2C_LOY_TFx"],
    ["TFx-adjusted LOY, median (range)", f"{pos.TFx_adjusted_LOY_primary_pct.median():.1f}% ({pos.TFx_adjusted_LOY_primary_pct.min():.1f}–{pos.TFx_adjusted_LOY_primary_pct.max():.1f}%)", "", "Fig2C_LOY_TFx"],
    ["Lowest alternative TFx-adjusted estimate", f"{pos.TFx_adjusted_min_pct.min():.1f}%", "", "Fig2C_LOY_TFx"],
    ["Matched blood LOY, 9 evaluable (sWGS)", f"{pr.matched_blood_LOY_pct.min():.1f} to {pr.matched_blood_LOY_pct.max():.1f}%", "", "Fig3A_tumor_vs_blood"],
    ["Median paired difference (7 LOY+)", f"{(pp.tumor_LOY_pct - pp.matched_blood_LOY_pct).median():.1f} points", f"exact Wilcoxon P = {w.pvalue:.3f}", "Fig3A_tumor_vs_blood"],
    ["Blood baseline chrY:chrX (10 non-LOY)", "0.952", "SD 0.0205; CV 2.2%", "Fig3B_parameters"],
    ["USC-027 blood", "73.9% (sWGS); 71.6% (dPCR)", "age 97; ~34 SD below baseline (descriptive)", "Fig3B_calibration"],
    ["USC-027 / USC-051 tumor chrY:chrX", f"{s1.chrY_chrX[0]:.2f} / {s1.chrY_chrX[1]:.2f}", "no relative Y loss", "FigS1_summary"],
    ["dPCR cohort", "61 men; median age 68 (21–97)", "", "Fig4AB_dPCR_cohort"],
    ["dPCR mLOY median (IQR)", f"{co.mLOY_pct.median():.1f}% ({co.mLOY_pct.quantile(.25):.1f} to {co.mLOY_pct.quantile(.75):.1f}%)", "", "Fig4AB_dPCR_cohort"],
    ["dPCR above 13.2% threshold", "4/61", f"6.6% (95% CI {k2[0]:.1f}–{k2[1]:.1f}%)", "Fig4AB_dPCR_cohort"],
    ["dPCR median CI half-width", f"{hw:.1f} points", "", "Fig4AB_dPCR_cohort"],
    ["mLOY vs age", f"rho = {rho.statistic:.2f}, P = {rho.pvalue:.3f}", f"excluding 027: rho = {rho2.statistic:.2f}, P = {rho2.pvalue:.3f}", "Fig4AB_dPCR_cohort"],
    ["dPCR mLOY, 9 evaluable patients", f"{ev.blood_mLOY_dPCR_pct.min():.1f} to {ev.blood_mLOY_dPCR_pct.max():.1f}%", "all below 13.2%", "Fig4C_concordance"],
    ["sWGS–dPCR concordance", f"{int(cc.concordant.sum())}/{len(cc)}", "", "Fig4C_concordance"],
    ["Bland–Altman bias (LoA)", f"{sa.difference_sWGS_minus_dPCR.mean():.1f} ({sa.difference_sWGS_minus_dPCR.mean()-1.96*sa.difference_sWGS_minus_dPCR.std():.1f} to {sa.difference_sWGS_minus_dPCR.mean()+1.96*sa.difference_sWGS_minus_dPCR.std():.1f})", "percentage points", "FigS3_agreement"],
]
key = pd.DataFrame(K, columns=["Quantity", "Value", "Detail", "Source sheet"])
print(key.to_string())

order = ["README", "Index", "Key_numbers", "Fig1_cohort"] + [k for k in x if k.startswith("Fig2") or k.startswith("Fig3") or k.startswith("Fig4")] \
        + ["FigS1_summary", "FigS1_seg_027_051"] + [k for k in x if k.startswith("FigS")] + ["dPCR_female_controls", "blood_sWGS_summary", "tumor_sWGS_summary"]
data = dict(x); data.update(README=r, Index=idx, Key_numbers=key, Fig1_cohort=fig1, FigS1_summary=s1.round(4), FigS1_seg_027_051=seg)
with pd.ExcelWriter(F, engine="openpyxl") as wr:
    for k in order: data[k].to_excel(wr, sheet_name=k, index=False)
wb = load_workbook(F)
for ws in wb.worksheets:
    for c in ws[1]: c.font = Font(bold=True, color="FFFFFF"); c.fill = PatternFill("solid", fgColor="4A5A6A"); c.alignment = Alignment(wrap_text=True, vertical="center")
    ws.freeze_panes = "A2"
    for col in ws.columns:
        L = max(len(str(c.value)) if c.value is not None else 0 for c in list(col)[:200])
        ws.column_dimensions[col[0].column_letter].width = min(max(10, L + 2), 90)
wb.save(F); print(wb.sheetnames)
