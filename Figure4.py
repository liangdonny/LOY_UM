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
CUT = 0.10
EXCL = ["027", "051"]

def yx(df, cols):
    y = df[(df.chrom == 24) & (df.chrompos >= 2.8e6) & (df.chrompos <= 26.6e6)][cols].mean()
    return y / df[df.chrom == 23][cols].mean()

b = pd.read_csv(DATA / "blood.csv"); bc = [c for c in b.columns if c.startswith("USC_")]
B = yx(b, bc); B.index = [c.split("_")[2] for c in B.index]
raw = pd.read_excel(DATA / "seg.xlsx", "raw_seg"); tc = [c for c in raw.columns if c.startswith("UM_")]
T = yx(raw, tc); T.index = [c.split("_")[1] for c in T.index]

null = B.drop("027"); n, BASE, sd = len(null), null.mean(), null.std(ddof=1)
bl = (1 - B / BASE) * 100; tu = (1 - T / BASE) * 100
ktol = lambda p: stats.nct.ppf(0.95, n - 1, stats.norm.ppf(p) * np.sqrt(n)) / np.sqrt(n)
tl99, tl999 = ktol(0.99) * sd / BASE * 100, ktol(0.999) * sd / BASE * 100

ev = tu.index.difference(EXCL)
e = pd.DataFrame({"tumor": tu[ev], "blood": bl[ev]}).sort_values("tumor", ascending=False)
e["pos"] = e.tumor >= CUT * 100

fig, (a, bx, c) = plt.subplots(1, 3, figsize=(7.4, 3.0), gridspec_kw=dict(width_ratios=[1.05, 0.85, 1.05], wspace=0.55))

# A: paired tumor vs matched blood
for s, r in e.iterrows():
    col = C_LOY if r.pos else C_NEG
    a.plot([0, 1], [r.tumor, r.blood], color=col, lw=0.9, alpha=0.8)
    a.scatter(0, r.tumor, s=22, color=col, zorder=3)
    a.scatter(1, r.blood, s=22, color=C_BLOOD, zorder=3)
labs = []
for s, r in e.iterrows():
    if labs and abs(labs[-1][0] - r.tumor) < 3: labs[-1][1].append(s)
    else: labs.append([r.tumor, [s]])
for v, nm in labs:
    a.text(-0.07, v, ", ".join(sorted(nm)), ha="right", va="center", fontsize=5.8, color="0.35")
a.axhline(CUT * 100, color="0.5", lw=0.7, ls="--")
a.set_xlim(-0.65, 1.25); a.set_ylim(-8, 100)
a.set_xticks([0, 1]); a.set_xticklabels(["Tumor", "Matched\nblood"])
a.set_ylabel("LOY fraction (%)")
a.legend(handles=[Line2D([], [], marker="o", color=C_LOY, ms=4, lw=0.9, label="Tumor LOY+"),
                  Line2D([], [], marker="o", color=C_NEG, ms=4, lw=0.9, label="Tumor LOY−"),
                  Line2D([], [], marker="o", ls="", color=C_BLOOD, ms=4, label="Blood")],
         frameon=False, fontsize=6.2, loc="upper right", bbox_to_anchor=(1.12, 1.0), handlelength=1.4)
a.set_title("A", loc="left", fontsize=9, fontweight="bold", x=-0.12)

# B: assay calibration — non-LOY blood noise and positive control
rng = np.random.default_rng(3)
bx.scatter(rng.uniform(-0.15, 0.15, n), bl.drop("027"), s=18, color=C_BLOOD, zorder=3)
bx.scatter([1], [bl["027"]], s=24, color=C_LOY, zorder=3)
bx.text(1, bl["027"] - 7, "027\ndPCR 72%", ha="center", va="top", fontsize=6, color="0.35")
bx.axhspan(tl99, tl999, color="0.9", zorder=0)
bx.axhline(CUT * 100, color="0.5", lw=0.7, ls="--")
bx.text(1.55, (tl99 + tl999) / 2 + 4.5, "99–99.9%\ntolerance\nlimits", ha="right", va="bottom", fontsize=5.6, color="0.45")
bx.set_xlim(-0.5, 1.6); bx.set_ylim(-8, 100)
bx.set_xticks([0, 1]); bx.set_xticklabels([f"Non-LOY\nblood\n(n = {n})", "Positive\ncontrol"])
bx.set_ylabel("Blood LOY fraction (%)")
bx.set_title("B", loc="left", fontsize=9, fontweight="bold", x=-0.12)

# C: threshold invariance
ts = np.linspace(0, 30, 301)
c.step(ts, [(e.tumor >= t).sum() for t in ts], where="post", color=C_LOY, lw=1.2, label=f"Tumors (of {len(e)})")
c.step(ts, [(bl.drop("027") >= t).sum() for t in ts], where="post", color=C_BLOOD, lw=1.2, label=f"Blood (of {n})")
gap_lo = max(e.tumor[~e.pos].max(), bl.drop("027").max()); gap_hi = e.tumor[e.pos].min()
c.axvspan(gap_lo, gap_hi, color="0.92", zorder=0)
c.text((gap_lo + gap_hi) / 2 + 2, 9.6, f"identical calls\n{gap_lo:.1f}–{gap_hi:.1f}%", ha="center", va="top", fontsize=6.2, color="0.35")
c.axvline(CUT * 100, color="0.5", lw=0.7, ls="--")
c.set_xlim(0, 30); c.set_ylim(0, 10.5); c.set_yticks(range(0, 11, 2))
c.set_xlabel("LOY threshold (%)"); c.set_ylabel("Samples called LOY-positive")
c.legend(frameon=False, fontsize=6.2, loc="center right", bbox_to_anchor=(1.0, 0.45))
c.set_title("C", loc="left", fontsize=9, fontweight="bold", x=-0.12)

fig.savefig(FIGS / "main" / "Figure4.png", dpi=600, bbox_inches="tight")
fig.savefig(FIGS / "main" / "Figure4.pdf", bbox_inches="tight")
print(f"BASE={BASE:.4f} SD={sd:.4f} TL99={tl99:.1f} TL99.9={tl999:.1f} gap={gap_lo:.1f}-{gap_hi:.1f}")
print(e.round(1).to_string()); print("027 blood", round(bl['027'], 1))
