import pandas as pd, numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from pathlib import Path
ROOT = Path("/Users/donny/LOY_UM")  # local clone of this repository; edit if cloned elsewhere
DATA, RESULTS, FIGS = ROOT / "data", ROOT / "results", ROOT / "figures"
for _d in (RESULTS / "Source_Data", FIGS / "main", FIGS / "supplementary"): _d.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6})
C_LOY, C_NEG, C_BLOOD = "#B5545C", "#9A9A9A", "#6F8BA4"
PAIRED = [5, 27, 33, 45, 47, 48, 51, 64, 68, 71, 140]
EVAL = [5, 33, 45, 47, 48, 64, 68, 71, 140]

m = pd.read_csv(DATA / "dpcr_clean.csv", index_col=0)
over = m.ratio[m.ratio > 1]; sd = np.sqrt(((over - 1) ** 2).mean()); CUTR = 1 - 2 * sd; CUT = (1 - CUTR) * 100
m["pos"] = m.ratio < CUTR

# sWGS (blood and tumor), baseline = mean of 10 non-LOY matched bloods
def yx(df, cols):
    y = df[(df.chrom == 24) & (df.chrompos >= 2.8e6) & (df.chrompos <= 26.6e6)][cols].mean()
    return y / df[df.chrom == 23][cols].mean()
b = pd.read_csv(DATA / "blood.csv"); bc = [c for c in b.columns if c.startswith("USC_")]
B = yx(b, bc); B.index = [int(c.split("_")[2]) for c in B.index]
BASE = B.drop(27).mean(); sw_b = (1 - B / BASE) * 100
raw = pd.read_excel(DATA / "seg.xlsx", "raw_seg"); tc = [c for c in raw.columns if c.startswith("UM_")]
T = yx(raw, tc); T.index = [int(c.split("_")[1]) for c in T.index]; sw_t = (1 - T / BASE) * 100

fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.0), gridspec_kw=dict(hspace=0.5, wspace=0.32, width_ratios=[1.35, 1]))
a, bx, c, dx = axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]

# A: ranked mLOY with 95% CI
o = m.sort_values("mloy", ascending=False).reset_index()
x = np.arange(len(o))
col = [C_LOY if p else C_BLOOD for p in o.pos]
a.vlines(x, o.lo, o.hi, color=col, lw=0.6, alpha=0.6)
a.scatter(x, o.mloy, s=11, color=col, zorder=3)
pi = o.index[o.id.isin(PAIRED)]
a.scatter(pi, o.mloy[pi], s=34, facecolor="none", edgecolor="0.15", lw=0.8, zorder=4)
a.axhline(CUT, color="0.5", lw=0.7, ls="--"); a.axhline(0, color="0.85", lw=0.6)
i27 = o.index[o.id == 27][0]
a.annotate("027", (i27, o.mloy[i27]), xytext=(i27 + 3, o.mloy[i27] - 2), fontsize=6.5, color="0.3", va="center")
a.set_xlim(-1.5, len(o) + 0.5); a.set_ylim(-28, 82); a.set_xticks([])
a.set_xlabel(f"Male patients ranked by mLOY (n = {len(o)})"); a.set_ylabel("Blood mLOY fraction (%)")
a.legend(handles=[Line2D([], [], marker="o", ls="", color=C_LOY, ms=4, label=f"mLOY+ ({int(m.pos.sum())}/{len(m)})"),
                  Line2D([], [], marker="o", ls="", color=C_BLOOD, ms=4, label="mLOY−"),
                  Line2D([], [], marker="o", ls="", mfc="none", mec="0.15", ms=6, label="Patients with tumor sWGS"),
                  Line2D([], [], color="0.5", lw=0.7, ls="--", label=f"Threshold ({CUT:.1f}%)")],
         frameon=False, fontsize=6.2, loc="upper right", handlelength=1.4)
a.set_title("A", loc="left", fontsize=9, fontweight="bold", x=-0.12)

# B: mLOY vs age
rho, p = stats.spearmanr(m.age, m.mloy)
bx.scatter(m.age, m.mloy, s=12, color=[C_LOY if q else C_BLOOD for q in m.pos], zorder=3)
bx.axhline(CUT, color="0.5", lw=0.7, ls="--"); bx.axhline(0, color="0.85", lw=0.6)
bx.set_xlabel("Age at diagnosis (years)"); bx.set_ylabel("Blood mLOY fraction (%)")
bx.set_ylim(-28, 82); bx.set_xlim(15, 100)
bx.text(0.03, 0.97, f"Spearman ρ = {rho:.2f}, P = {p:.3f}", transform=bx.transAxes, fontsize=6.5, va="top", color="0.3")
bx.set_title("B", loc="left", fontsize=9, fontweight="bold", x=-0.12)

# C: cross-platform concordance in matched blood
cc = pd.DataFrame({"swgs": sw_b.reindex(PAIRED), "dpcr": m.mloy.reindex(PAIRED)})
c.axvline(10, color="0.5", lw=0.7, ls="--"); c.axhline(CUT, color="0.5", lw=0.7, ls=":")
c.plot([-15, 85], [-15, 85], color="0.8", lw=0.6)
c.scatter(cc.swgs, cc.dpcr, s=22, color=[C_LOY if s == 27 else C_BLOOD for s in cc.index], zorder=3)
c.annotate("027", (cc.swgs[27], cc.dpcr[27]), xytext=(cc.swgs[27] - 3, cc.dpcr[27]), fontsize=6.5, color="0.3", ha="right", va="center")
both_pos = ((cc.swgs >= 10) & (cc.dpcr > CUT)).sum(); both_neg = ((cc.swgs < 10) & (cc.dpcr <= CUT)).sum()
c.text(0.97, 0.05, f"Concordant calls: {both_pos + both_neg}/{len(cc)}", transform=c.transAxes, ha="right", fontsize=6.5, color="0.3")
c.set_xlim(-15, 85); c.set_ylim(-15, 85)
c.set_xlabel("Blood LOY fraction, sWGS (%)"); c.set_ylabel("Blood mLOY fraction, dPCR (%)")
c.legend(handles=[Line2D([], [], color="0.5", lw=0.7, ls="--", label="sWGS threshold (10%)"),
                  Line2D([], [], color="0.5", lw=0.7, ls=":", label=f"dPCR threshold ({CUT:.1f}%)")],
         frameon=False, fontsize=6.2, loc="upper left", handlelength=1.8)
c.set_title("C", loc="left", fontsize=9, fontweight="bold", x=-0.12)

# D: tumor LOY (raw and TFx-adjusted) vs blood (summary)
rng = np.random.default_rng(2)
adj = pd.read_excel(RESULTS / "Source_Data" / "LOY_UM_MasterData.xlsx", "Fig2C_LOY_TFx")
adj.index = [int(p.split("-")[1]) for p in adj.Patient]
traw = sw_t.reindex(EVAL); tadj = adj.TFx_adjusted_LOY_primary_pct.reindex(EVAL)
jt = rng.uniform(-0.12, 0.12, len(EVAL))
for k, pid in enumerate(EVAL):
    pos = traw.loc[pid] >= 10; c = C_LOY if pos else C_NEG
    dx.plot([0 + jt[k], 1 + jt[k]], [traw.loc[pid], tadj.loc[pid]], color=c, lw=0.6, alpha=0.5, zorder=2)
    dx.scatter(0 + jt[k], traw.loc[pid], s=14, facecolor="white", edgecolor=c, lw=1, zorder=3)
    dx.scatter(1 + jt[k], tadj.loc[pid], s=14, color=c, zorder=3)
dx.hlines(np.median(traw), -0.28, 0.28, color="0.2", lw=1); dx.hlines(np.median(tadj), 0.72, 1.28, color="0.2", lw=1)
groups = [("Tumor\nraw", traw, C_LOY), ("Tumor\nTFx-adj.", tadj, C_LOY), ("Matched\nblood", sw_b.reindex(EVAL), C_BLOOD),
          ("Blood\ncohort", m.mloy, C_BLOOD)]
for i, (lab, v, colr) in list(enumerate(groups))[2:]:
    jit = rng.uniform(-0.18, 0.18, len(v))
    cols = [C_LOY if (i == 3 and val > CUT) else colr for val in v]
    dx.scatter(i + jit, v, s=11, color=cols, zorder=3, alpha=0.9)
    dx.hlines(np.median(v), i - 0.28, i + 0.28, color="0.2", lw=1)
dx.axvline(1.5, color="0.85", lw=0.6)
dx.text(0.5, 103, "sWGS", ha="center", fontsize=6.3, color="0.4"); dx.text(2, 103, "sWGS", ha="center", fontsize=6.3, color="0.4"); dx.text(3, 103, "dPCR", ha="center", fontsize=6.3, color="0.4")
dx.axhline(0, color="0.85", lw=0.6)
dx.set_xticks(range(4)); dx.set_xticklabels([g[0] for g in groups], fontsize=6.5)
dx.set_ylim(-28, 108); dx.set_ylabel("LOY fraction (%)")
dx.set_title("D", loc="left", fontsize=9, fontweight="bold", x=-0.12)

fig.savefig(FIGS / "main" / "Figure5.png", dpi=600, bbox_inches="tight")
fig.savefig(FIGS / "main" / "Figure5.pdf", bbox_inches="tight")
print(f"SD={sd:.4f} cut ratio={CUTR:.4f} mLOY cut={CUT:.1f}% pos={m.pos.sum()}/{len(m)}")
print(f"rho={rho:.3f} p={p:.4f}"); print(cc.round(1)); print("medians", {g[0]: round(float(np.median(g[1])), 1) for g in groups})
