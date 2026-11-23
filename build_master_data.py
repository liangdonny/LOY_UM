"""Build LOY_UM_MasterData.xlsx: one sheet per figure panel, all derived from the raw inputs."""
import pandas as pd, numpy as np
from scipy import stats

from pathlib import Path
ROOT = Path("/Users/donny/LOY_UM")  # local clone of this repository; edit if cloned elsewhere
DATA, RESULTS, FIGS = ROOT / "data", ROOT / "results", ROOT / "figures"
for _d in (RESULTS / "Source_Data", FIGS / "main", FIGS / "supplementary"): _d.mkdir(parents=True, exist_ok=True)
SEG = DATA / "seg.xlsx"            # tumor segments (TFx-corrected, raw, TFx sheets)
BLOOD = DATA / "blood.csv"         # matched blood segments (11)
PCR = DATA / "LOY_PCR.csv"         # QIAcuity export
HIST = DATA / "UM_karyotype_pooled_Horsman1993_Prescher1990.csv"  # written by FigureS4.py
OUT = RESULTS / "Source_Data" / "LOY_UM_MasterData.xlsx"
CUT = 10.0
pid = lambda s: "USC-" + str(s).zfill(3)

def yx_parts(df, cols):
    win = (df.chrom == 24) & (df.chrompos >= 2.8e6) & (df.chrompos <= 26.6e6)
    return df.loc[win, cols].mean(), df.loc[df.chrom == 23, cols].mean(), df.loc[df.chrom <= 22, cols].mean()

# ---------- blood sWGS ----------
b = pd.read_csv(BLOOD); bc = [c for c in b.columns if c.startswith("USC_")]
bY, bX, bA = yx_parts(b, bc)
B = pd.DataFrame({"chrY_mean": bY, "chrX_mean": bX, "autosome_mean": bA}); B.index = [c.split("_")[2] for c in B.index]
B["chrY:chrX"] = B.chrY_mean / B.chrX_mean; B["2chrY:autosome"] = 2 * B.chrY_mean / B.autosome_mean
null = B.drop("027"); BASE, SD = null["chrY:chrX"].mean(), null["chrY:chrX"].std(ddof=1); BASE_A = null["2chrY:autosome"].mean()
B["LOY_fraction_pct"] = (1 - B["chrY:chrX"] / BASE) * 100
B["LOY_fraction_YA_pct"] = (1 - B["2chrY:autosome"] / BASE_A) * 100

# ---------- tumor sWGS ----------
raw = pd.read_excel(SEG, "raw_seg"); tc = [c for c in raw.columns if c.startswith("UM_")]
tY, tX, _ = yx_parts(raw, tc)
st = pd.read_excel(SEG, "Supp_TF_sensitivity", header=None); h = st.index[st[0] == "Sample"][0]
tt = st.iloc[h + 1:h + 1 + len(tc)].copy(); tt.columns = st.iloc[h].tolist(); tt.index = [s.split("_")[1] for s in tt.Sample]
T = pd.DataFrame({"chrY_mean": tY.values, "chrX_mean": tX.values}, index=[c.split("_")[1] for c in tc])
T["chrY:chrX"] = T.chrY_mean / T.chrX_mean
T["raw_LOY_fraction_pct"] = (1 - T["chrY:chrX"] / BASE) * 100
T["LOY_call"] = np.where(T.raw_LOY_fraction_pct >= CUT, "LOY+", "LOY-")
for col, lab in [("TF – primary", "TFx_primary"), ("TF – alt 1", "TFx_alt1"), ("TF – alt 2", "TFx_alt2")]:
    T[lab] = pd.to_numeric(tt[col], errors="coerce").reindex(T.index)
for lab in ["primary", "alt1", "alt2"]:
    T[f"TFx_adjusted_LOY_{lab}_pct"] = (T.raw_LOY_fraction_pct / T[f"TFx_{lab}"]).clip(0, 100)
adj = T[[c for c in T.columns if c.startswith("TFx_adjusted")]]
T["TFx_adjusted_min_pct"], T["TFx_adjusted_max_pct"] = adj.min(axis=1), adj.max(axis=1)
T["primary_anchor"] = tt["Primary anchor"].reindex(T.index)
T["TFx_low_confidence"] = T.index.isin(["033", "048"])

# arm-level (TFx-corrected segments)
seg = pd.read_excel(SEG, "seg_TFcorrected_data_101526"); S = [c for c in seg.columns if c.startswith("UM_")]
cen = {1:125.0,3:91.0,6:61.0,8:45.6,16:36.6}
def arm(ch, a):
    m = (seg.chrom == ch) & ((seg.chrompos / 1e6 < cen[ch]) if a == "p" else (seg.chrompos / 1e6 >= cen[ch]))
    return seg.loc[m, S].mean()
ARM = pd.DataFrame({k: arm(int(k[:-1]), k[-1]) for k in ["1p", "3p", "3q", "6p", "6q", "8p", "8q", "16q"]}); ARM.index = [s.split("_")[1] for s in S]
order = list(T[T.LOY_call == "LOY+"].sort_values("raw_LOY_fraction_pct", ascending=False).index) + list(T[T.LOY_call == "LOY-"].sort_values("raw_LOY_fraction_pct", ascending=False).index)
ev = {"LOY": T.LOY_call == "LOY+", "-3": (ARM["3p"] < 0.75) | (ARM["3q"] < 0.75), "+8q": ARM["8q"] > 1.25, "+6p": ARM["6p"] > 1.25,
      "-8p": ARM["8p"] < 0.75, "-1p": ARM["1p"] < 0.75, "-16q": ARM["16q"] < 0.75, "-6q": ARM["6q"] < 0.75}
EV = pd.DataFrame(ev)

# genome-wide TFx-adjusted copy state (long)
normal = np.where(seg.chrom >= 23, 0.5, 1.0)
gw = seg[["chrom", "chrompos"]].copy()
for s in S: gw[pid(s.split("_")[1])] = seg[s] / normal
gw = gw[~((gw.chrom == 24) & ((gw.chrompos < 2.8e6) | (gw.chrompos > 26.6e6)))]

# ---------- dPCR ----------
d = pd.read_csv(PCR, encoding="latin1")
d.columns = ["id", "sex", "age", "scna", "loypos", "mix", "target", "conc", "ci", "pv", "pp", "pn", "thr", "vol", "ref", "ratio"]
d = d[d.id != "NTC"]
y = d[d.target == "AMELY"].set_index("id"); x = d[d.target == "AMELX"].set_index("id")
P = pd.DataFrame({"sex": y.sex, "age_at_diagnosis": y.age.astype(float), "AMELY_cp_per_uL": y.conc.astype(float), "AMELX_cp_per_uL": x.conc.astype(float),
                  "AMELY_rel_CI95": pd.to_numeric(y.ci, errors="coerce"), "AMELX_rel_CI95": pd.to_numeric(x.ci, errors="coerce"),
                  "AMELY_pos_partitions": y.pp.astype(int), "AMELX_pos_partitions": x.pp.astype(int), "valid_partitions": y.pv.astype(int)})
P.index = [str(i).zfill(3) for i in P.index]
FEM = P[P.sex == "F"].copy(); P = P[P.sex == "M"].copy()
P["AMELY:AMELX"] = P.AMELY_cp_per_uL / P.AMELX_cp_per_uL
rel = np.sqrt(P.AMELY_rel_CI95 ** 2 + P.AMELX_rel_CI95 ** 2)
P["mLOY_pct"] = (1 - P["AMELY:AMELX"]) * 100
P["mLOY_CI95_low_pct"] = (1 - P["AMELY:AMELX"] * (1 + rel)) * 100; P["mLOY_CI95_high_pct"] = (1 - P["AMELY:AMELX"] * (1 - rel)) * 100
over = P["AMELY:AMELX"][P["AMELY:AMELX"] > 1]; DSD = np.sqrt(((over - 1) ** 2).mean()); DCUT = 1 - 2 * DSD
P["mLOY_call"] = np.where(P["AMELY:AMELX"] < DCUT, "mLOY+", "mLOY-")
PAIRED = ["005", "027", "033", "045", "047", "048", "051", "064", "068", "071", "140"]
P["has_tumor_sWGS"] = P.index.isin(PAIRED)
P["flag"] = np.where(P.index == "046", "AMELY and AMELX values identical; verify raw export", "")
rho, pval = stats.spearmanr(P.age_at_diagnosis, P.mLOY_pct)

EVAL = [s for s in T.index if s not in ("027", "051")]

# ---------- sheets ----------
sheets = {}
readme = pd.DataFrame({"Item": [
    "Project", "Analytic tumor cohort", "Excluded tumors", "sWGS chrY window", "sWGS chrY:chrX", "sWGS baseline",
    "sWGS LOY fraction", "sWGS LOY threshold", "TFx", "TFx-adjusted LOY fraction", "dPCR mLOY", "dPCR threshold",
    "Arm calls", "Sample IDs", "Inputs", "Data flag"],
  "Definition": [
    "Tumor and host loss of chromosome Y in uveal melanoma (Paper 1)",
    "9 male intraocular fluid specimens (one per patient): " + ", ".join(pid(s) for s in EVAL),
    "USC-027, USC-051: ploidy-indeterminate complex hyperdiploid profiles (blood retained; 027 blood = positive control)",
    "Male-specific region, chrY 2.8–26.6 Mb (hg19)",
    "Mean segmented copy ratio over the chrY window ÷ mean over chrX",
    f"{BASE:.4f} = mean chrY:chrX of 10 matched blood samples without LOY (all matched bloods except USC-027); SD {SD:.4f}",
    "1 − chrY:chrX ÷ baseline, in %",
    "≥10% (one-sided normal tolerance limit covering 99.7% of non-LOY bloods, 95% confidence)",
    "Tumor fraction from one clonal autosomal anchor: chr3 monosomy, else clean single-copy loss, else gain",
    "Raw LOY fraction ÷ TFx, capped at 100%",
    "1 − AMELY:AMELX (theoretical baseline 1.0); 95% CI by delta method from channel CIs",
    f"AMELY:AMELX < {DCUT:.4f} (2 SD below 1.0; SD {DSD:.4f} = RMS of ratios > 1.0, n = {len(over)}) = mLOY > {(1 - DCUT) * 100:.1f}%",
    "TFx-adjusted arm mean < 0.75 = loss, > 1.25 = gain; hg19 centromeres",
    "USC-### = patient; one specimen per patient",
    "seg_TFcorrected_data_101526.xlsx; seg_data_blood_11pairs.csv; LOY_PCR.csv; Horsman 1993 and Prescher 1990 tables",
    "dPCR sample 046: AMELY and AMELX concentrations identical"]})
sheets["README"] = readme

sheets["Fig2A_frequency"] = pd.DataFrame({"alteration": EV.columns, "n_tumors": EV.loc[EVAL].sum().values, "n_evaluable": len(EVAL)})
f2b = gw[["chrom", "chrompos"] + [pid(s) for s in order]]
sheets["Fig2B_genomewide"] = f2b
c2 = T.loc[order, ["chrY:chrX", "raw_LOY_fraction_pct", "LOY_call", "TFx_primary", "TFx_alt1", "TFx_alt2", "TFx_adjusted_LOY_primary_pct",
                   "TFx_adjusted_min_pct", "TFx_adjusted_max_pct", "primary_anchor", "TFx_low_confidence"]].copy()
c2.index = [pid(s) for s in c2.index]; c2.index.name = "Patient"; sheets["Fig2C_LOY_TFx"] = c2.reset_index()
a2 = ARM.loc[order].copy(); a2.index = [pid(s) for s in a2.index]; a2.index.name = "Patient"; sheets["Fig2_arm_copy_state"] = a2.reset_index()

f3a = pd.DataFrame({"Patient": [pid(s) for s in order], "tumor_LOY_pct": T.loc[order, "raw_LOY_fraction_pct"].values,
                    "tumor_call": T.loc[order, "LOY_call"].values, "matched_blood_LOY_pct": B.loc[order, "LOY_fraction_pct"].values})
sheets["Fig3A_tumor_vs_blood"] = f3a
pos = f3a[f3a.tumor_call == "LOY+"]; w = stats.wilcoxon(pos.tumor_LOY_pct, pos.matched_blood_LOY_pct)
k = lambda p: stats.nct.ppf(0.95, len(null) - 1, stats.norm.ppf(p) * np.sqrt(len(null))) / np.sqrt(len(null))
f3b = B[["chrY:chrX", "LOY_fraction_pct"]].copy(); f3b["role"] = np.where(f3b.index == "027", "positive control", "non-LOY reference")
f3b.index = [pid(s) for s in f3b.index]; f3b.index.name = "Patient"
sheets["Fig3B_calibration"] = f3b.reset_index()
sheets["Fig3B_parameters"] = pd.DataFrame({"parameter": ["n non-LOY bloods", "mean chrY:chrX", "SD", "CV %", "mean - 2SD LOY %", "mean - 3SD LOY %",
                                                          "tolerance 99%/95% LOY %", "tolerance 99.9%/95% LOY %", "Wilcoxon 7 LOY+ pairs P", "median paired difference (points)"],
    "value": [len(null), BASE, SD, SD / BASE * 100, 2 * SD / BASE * 100, 3 * SD / BASE * 100, k(0.99) * SD / BASE * 100, k(0.999) * SD / BASE * 100,
              w.pvalue, (pos.tumor_LOY_pct - pos.matched_blood_LOY_pct).median()]})
ts = np.round(np.arange(0, 30.1, 0.5), 1)
sheets["Fig3C_threshold_sweep"] = pd.DataFrame({"threshold_pct": ts,
    "tumors_LOY+_of_9": [(T.loc[EVAL, "raw_LOY_fraction_pct"] >= t).sum() for t in ts],
    "nonLOY_bloods_LOY+_of_10": [(null.index.map(lambda s: B.loc[s, "LOY_fraction_pct"]) >= t).sum() for t in ts]})

f4 = P[["age_at_diagnosis", "AMELY_cp_per_uL", "AMELX_cp_per_uL", "AMELY_rel_CI95", "AMELX_rel_CI95", "AMELY:AMELX", "mLOY_pct",
        "mLOY_CI95_low_pct", "mLOY_CI95_high_pct", "mLOY_call", "has_tumor_sWGS", "flag"]].copy()
f4.index = [pid(s) for s in f4.index]; f4.index.name = "Patient"; sheets["Fig4AB_dPCR_cohort"] = f4.reset_index().sort_values("mLOY_pct", ascending=False)
sheets["Fig4B_stats"] = pd.DataFrame({"parameter": ["n men", "median age", "age range", "Spearman rho (age vs mLOY)", "P", "mLOY+ n", "median mLOY %"],
    "value": [len(P), P.age_at_diagnosis.median(), f"{P.age_at_diagnosis.min():.0f}-{P.age_at_diagnosis.max():.0f}", rho, pval, int((P.mLOY_call == "mLOY+").sum()), P.mLOY_pct.median()]})
f4c = pd.DataFrame({"Patient": [pid(s) for s in PAIRED], "blood_LOY_sWGS_pct": B.loc[PAIRED, "LOY_fraction_pct"].values,
                    "blood_mLOY_dPCR_pct": P.loc[PAIRED, "mLOY_pct"].values,
                    "call_sWGS": np.where(B.loc[PAIRED, "LOY_fraction_pct"] >= CUT, "LOY+", "LOY-"), "call_dPCR": P.loc[PAIRED, "mLOY_call"].values})
f4c["concordant"] = (f4c.call_sWGS == "LOY+") == (f4c.call_dPCR == "mLOY+"); sheets["Fig4C_concordance"] = f4c
sheets["Fig4D_groups"] = pd.concat([
    pd.DataFrame({"group": "Tumor (sWGS)", "Patient": [pid(s) for s in EVAL], "LOY_pct": T.loc[EVAL, "raw_LOY_fraction_pct"].values}),
    pd.DataFrame({"group": "Matched blood (sWGS)", "Patient": [pid(s) for s in EVAL], "LOY_pct": B.loc[EVAL, "LOY_fraction_pct"].values}),
    pd.DataFrame({"group": "Blood cohort (dPCR)", "Patient": [pid(s) for s in P.index], "LOY_pct": P.mLOY_pct.values})])
sheets["dPCR_female_controls"] = FEM[["AMELY_cp_per_uL", "AMELX_cp_per_uL"]].rename_axis("Sample").reset_index()

s2a = pd.DataFrame({"Patient": [pid(s) for s in EVAL], "LOY_pct_baseline_primary": B.loc[EVAL, "LOY_fraction_pct"].values,
                    "LOY_pct_baseline_0.936": (1 - B.loc[EVAL, "chrY:chrX"] / 0.9356).values * 100,
                    "LOY_pct_baseline_1.0": (1 - B.loc[EVAL, "chrY:chrX"]).values * 100})
sheets["FigS2A_blood_baselines"] = s2a
bgw = b[["chrom", "chrompos"]].copy()
for s in EVAL:
    col = f"USC_UM_{s}_{'Blood' if s == '005' else 'gDNA'}"; bgw[pid(s)] = b[col] / np.where(b.chrom >= 23, 0.5, 1.0)
sheets["FigS2B_blood_genomewide"] = bgw[~((bgw.chrom == 24) & ((bgw.chrompos < 2.8e6) | (bgw.chrompos > 26.6e6)))]
ba = f4c[["Patient", "blood_LOY_sWGS_pct", "blood_mLOY_dPCR_pct"]].copy()
ba["mean_pct"] = (ba.blood_LOY_sWGS_pct + ba.blood_mLOY_dPCR_pct) / 2; ba["difference_sWGS_minus_dPCR"] = ba.blood_LOY_sWGS_pct - ba.blood_mLOY_dPCR_pct
ba["LOY_sWGS_Yautosome_pct"] = B.loc[PAIRED, "LOY_fraction_YA_pct"].values
sheets["FigS3_agreement"] = ba
bias, sdd = ba.difference_sWGS_minus_dPCR.mean(), ba.difference_sWGS_minus_dPCR.std(ddof=1)
sheets["FigS3_parameters"] = pd.DataFrame({"parameter": ["bias (points)", "SD of differences", "lower 95% LoA", "upper 95% LoA", "2Y:autosome baseline"],
                                           "value": [bias, sdd, bias - 1.96 * sdd, bias + 1.96 * sdd, BASE_A]})
sheets["FigS4_historical"] = pd.read_csv(HIST)

braw = B[["chrY_mean", "chrX_mean", "autosome_mean", "chrY:chrX", "2chrY:autosome", "LOY_fraction_pct", "LOY_fraction_YA_pct"]].copy()
braw.index = [pid(s) for s in braw.index]; braw.index.name = "Patient"; sheets["blood_sWGS_summary"] = braw.reset_index()
traw = T[["chrY_mean", "chrX_mean", "chrY:chrX", "raw_LOY_fraction_pct"]].copy(); traw.index = [pid(s) for s in traw.index]; traw.index.name = "Patient"
sheets["tumor_sWGS_summary"] = traw.reset_index()

with pd.ExcelWriter(OUT) as wr:
    for name, df in sheets.items():
        df.to_excel(wr, sheet_name=name[:31], index=False)
print("written", OUT); print(f"BASE={BASE:.4f} SD={SD:.4f} DCUT={DCUT:.4f} rho={rho:.3f} P={pval:.4f} wilcoxon P={w.pvalue:.4f} bias={bias:.2f}")
print(sheets["Fig2A_frequency"].to_string(index=False))
