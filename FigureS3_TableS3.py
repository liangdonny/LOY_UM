import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from pathlib import Path
ROOT = Path("/Users/donny/LOY_UM")  # local clone of this repository; edit if cloned elsewhere
DATA, RESULTS, FIGS = ROOT / "data", ROOT / "results", ROOT / "figures"
for _d in (RESULTS / "Source_Data", FIGS / "main", FIGS / "supplementary"): _d.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6})
C_LOY, C_BLOOD = "#B5545C", "#6F8BA4"
IDS = [5, 27, 33, 45, 47, 48, 51, 64, 68, 71, 140]

b = pd.read_csv(DATA / "blood.csv"); bc = [c for c in b.columns if c.startswith("USC_")]
win = (b.chrom == 24) & (b.chrompos >= 2.8e6) & (b.chrompos <= 26.6e6)
Y = b.loc[win, bc].mean(); X = b.loc[b.chrom == 23, bc].mean(); A = b.loc[b.chrom <= 22, bc].mean()
t = pd.DataFrame({"chrY_mean": Y, "chrX_mean": X, "auto_mean": A})
t.index = [int(c.split("_")[2]) for c in t.index]
t["YX"] = t.chrY_mean / t.chrX_mean; t["2YA"] = 2 * t.chrY_mean / t.auto_mean
null = t.drop(27)
bYX, bYA = null.YX.mean(), null["2YA"].mean()
t["LOY_sWGS_YX"] = (1 - t.YX / bYX) * 100
t["LOY_sWGS_YA"] = (1 - t["2YA"] / bYA) * 100

m = pd.read_csv(DATA / "dpcr_clean.csv", index_col=0)
over = m.ratio[m.ratio > 1]; sd = np.sqrt(((over - 1) ** 2).mean()); CUTR = 1 - 2 * sd; CUT = (1 - CUTR) * 100
t = t.join(m[["ratio", "mloy", "lo", "hi"]].rename(columns={"ratio": "dPCR_ratio", "mloy": "mLOY_dPCR", "lo": "dPCR_CI_low", "hi": "dPCR_CI_high"}))
t["call_sWGS"] = np.where(t.LOY_sWGS_YX >= 10, "LOY+", "LOY−")
t["call_dPCR"] = np.where(t.dPCR_ratio < CUTR, "LOY+", "LOY−")
t["concordant"] = t.call_sWGS == t.call_dPCR
t = t.loc[IDS]

# Bland–Altman (sWGS − dPCR), LOY fraction in %
mean_ = (t.LOY_sWGS_YX + t.mLOY_dPCR) / 2; diff = t.LOY_sWGS_YX - t.mLOY_dPCR
bias, sdd = diff.mean(), diff.std(ddof=1); lo, hi = bias - 1.96 * sdd, bias + 1.96 * sdd

fig, (a, c) = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw=dict(wspace=0.38, width_ratios=[1.25, 1]))
a.axhline(0, color="0.85", lw=0.6)
a.axhline(bias, color="0.25", lw=1)
for v in (lo, hi): a.axhline(v, color="0.5", lw=0.7, ls="--")
cols = [C_LOY if s == 27 else C_BLOOD for s in t.index]
a.scatter(mean_, diff, s=22, color=cols, zorder=3)
for s in [27, 33]:
    a.annotate(f"{s:03d}", (mean_[s], diff[s]), xytext=(5, 4 if s == 27 else -9), textcoords="offset points", fontsize=6.5, color="0.3")
a.text(0.98, bias - 1.2, f"bias {bias:+.1f}", transform=a.get_yaxis_transform(), ha="right", va="top", fontsize=6.5, color="0.25")
a.text(0.98, hi + 0.8, f"upper 95% limit {hi:+.1f}", transform=a.get_yaxis_transform(), ha="right", va="bottom", fontsize=6.5, color="0.45")
a.text(0.98, lo - 0.8, f"lower 95% limit {lo:+.1f}", transform=a.get_yaxis_transform(), ha="right", va="top", fontsize=6.5, color="0.45")
a.set_xlim(-10, 80); a.set_ylim(-20, 20)
a.set_xlabel("Mean blood LOY fraction, sWGS and dPCR (%)"); a.set_ylabel("sWGS − dPCR (percentage points)")
a.set_title("A", loc="left", fontsize=9, fontweight="bold", x=-0.12)

# B: Y:X vs Y:autosome (sWGS), all 11 bloods
c.plot([-10, 80], [-10, 80], color="0.8", lw=0.6)
c.axvline(10, color="0.5", lw=0.7, ls="--"); c.axhline(10, color="0.5", lw=0.7, ls="--")
c.scatter(t.LOY_sWGS_YX, t.LOY_sWGS_YA, s=22, color=cols, zorder=3)
c.annotate("027", (t.LOY_sWGS_YX[27], t.LOY_sWGS_YA[27]), xytext=(-5, 0), textcoords="offset points", ha="right", va="center", fontsize=6.5, color="0.3")
c.set_xlim(-10, 80); c.set_ylim(-10, 80)
c.set_xlabel("Blood LOY, chrY:chrX reference (%)"); c.set_ylabel("Blood LOY, chrY:autosome reference (%)")
c.set_title("B", loc="left", fontsize=9, fontweight="bold", x=-0.12)
fig.savefig(FIGS / "supplementary" / "FigureS3.png", dpi=600, bbox_inches="tight")
fig.savefig(FIGS / "supplementary" / "FigureS3.pdf", bbox_inches="tight")

out = t[["chrY_mean", "chrX_mean", "auto_mean", "YX", "2YA", "LOY_sWGS_YX", "LOY_sWGS_YA", "dPCR_ratio", "mLOY_dPCR", "dPCR_CI_low", "dPCR_CI_high", "call_sWGS", "call_dPCR", "concordant"]].copy()
out.index = [f"USC-{i:03d}" for i in out.index]; out.index.name = "Patient"
with pd.ExcelWriter(RESULTS / "TableS_crossplatform_blood.xlsx") as w:
    out.round(4).to_excel(w, sheet_name="per_patient")
    pd.DataFrame({"Parameter": ["sWGS chrY:chrX baseline (mean of 10 non-LOY bloods)", "sWGS 2Y:autosome baseline", "sWGS LOY threshold (%)", "dPCR SD (RMS of ratios >1.0)", "dPCR ratio threshold", "dPCR mLOY threshold (%)", "Bland–Altman bias (sWGS − dPCR, pp)", "95% limits of agreement (pp)", "Concordant calls"],
                  "Value": [round(bYX, 4), round(bYA, 4), 10, round(sd, 4), round(CUTR, 4), round(CUT, 1), round(bias, 1), f"{lo:.1f} to {hi:.1f}", f"{int(t.concordant.sum())}/{len(t)}"]}).to_excel(w, sheet_name="parameters", index=False)
print(out.round(2).to_string()); print(f"baseYX={bYX:.4f} base2YA={bYA:.4f} bias={bias:.2f} LoA {lo:.2f} {hi:.2f}")
print("max |YX-YA| excl 027:", (t.LOY_sWGS_YX - t.LOY_sWGS_YA).drop(27).abs().max().round(2))
