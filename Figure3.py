import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap

from pathlib import Path
ROOT = Path("/Users/donny/LOY_UM")  # local clone of this repository; edit if cloned elsewhere
DATA, RESULTS, FIGS = ROOT / "data", ROOT / "results", ROOT / "figures"
for _d in (RESULTS / "Source_Data", FIGS / "main", FIGS / "supplementary"): _d.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6})
C_LOY, C_LOSS, C_GAIN = "#B5545C", "#4A5A6A", "#E08A3C"
BASE, CUT = 0.9522, 0.10

seg = pd.read_excel(DATA / "seg.xlsx", "seg_TFcorrected_data_101526")
S = [c for c in seg.columns if c.startswith("UM_")]
st = pd.read_excel(DATA / "seg.xlsx", "Supp_TF_sensitivity", header=None)
h = st.index[st[0] == "Sample"][0]
L = st.iloc[h + 1:h + 1 + len(S)].copy(); L.columns = st.iloc[h].tolist(); L = L.set_index("Sample")
L = L[["Y:X ratio", "TF – primary", "TF – alt 1", "TF – alt 2"]].apply(pd.to_numeric, errors="coerce")
L.columns = ["yx", "tf", "tf1", "tf2"]
L["raw"] = 1 - L.yx / BASE
for a, b in [("adj", "tf"), ("adj1", "tf1"), ("adj2", "tf2")]:
    L[a] = (L.raw / L[b]).clip(0, 1)
L["pos"] = L.raw >= CUT
L["lo"] = L[["adj", "adj1", "adj2"]].min(axis=1); L["hi"] = L[["adj", "adj1", "adj2"]].max(axis=1)

cen = {1:125.0,2:93.3,3:91.0,4:50.4,5:48.4,6:61.0,7:59.9,8:45.6,9:49.0,10:40.2,11:53.7,12:35.8,13:17.9,14:17.6,15:19.0,16:36.6,17:24.0,18:17.2,19:26.5,20:27.5,21:13.2,22:14.7}
def arm(ch, a):
    m = (seg.chrom == ch) & ((seg.chrompos / 1e6 < cen[ch]) if a == "p" else (seg.chrompos / 1e6 >= cen[ch]))
    return seg.loc[m, S].mean()
A = pd.DataFrame({k: arm(int(k[:-1]), k[-1]) for k in ["1p", "3p", "3q", "6p", "6q", "8p", "8q", "16q"]})
ev = {"LOY": L.pos, "−3": (A["3p"] < 0.75) | (A["3q"] < 0.75), "+8q": A["8q"] > 1.25, "+6p": A["6p"] > 1.25,
      "−8p": A["8p"] < 0.75, "−1p": A["1p"] < 0.75, "−16q": A["16q"] < 0.75, "−6q": A["6q"] < 0.75}
F = pd.DataFrame(ev); n = len(S)
order = list(L[L.pos].sort_values("raw", ascending=False).index) + list(L[~L.pos].sort_values("raw", ascending=False).index)
lab = [s.replace("UM_", "") for s in order]; lowconf = {"UM_033", "UM_048"}; npos = int(L.pos.sum())

fig = plt.figure(figsize=(7.2, 3.3))
gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1], wspace=0.42)
# A
axA = fig.add_subplot(gs[0, 0]); fr = F.sum()
axA.bar(range(len(fr)), fr.values, width=0.62, color=[C_LOY if k == "LOY" else (C_GAIN if k.startswith("+") else C_LOSS) for k in fr.index])
for i, v in enumerate(fr.values): axA.text(i, v + 0.15, f"{v}/{n}", ha="center", va="bottom", fontsize=6)
axA.set_xticks(range(len(fr))); axA.set_xticklabels(fr.index, fontsize=6.3, rotation=0); axA.set_ylim(0, n + 0.8); axA.set_yticks([0, 3, 6, 9]); axA.set_ylabel("Tumors")
axA.legend(handles=[Patch(color=C_LOY, label="LOY"), Patch(color=C_LOSS, label="Arm loss"), Patch(color=C_GAIN, label="Arm gain")], frameon=False, fontsize=6.5, loc="upper right", handlelength=1)
axA.set_title("A", loc="left", fontsize=9, fontweight="bold", x=-0.16)


axC = fig.add_subplot(gs[0, 1])
# C dumbbell, same row order
for i, s in enumerate(order):
    y = i + 0.5; r = L.loc[s]
    if r.pos:
        axC.plot([r.lo * 100, r.hi * 100], [y, y], color=C_LOY, alpha=0.25, lw=5, solid_capstyle="round")
        axC.plot([r.raw * 100, r.adj * 100], [y, y], color=C_LOY, lw=1.1)
        axC.scatter(r.raw * 100, y, s=24, facecolor="white", edgecolor=C_LOY, lw=1.2, zorder=3)
        axC.scatter(r.adj * 100, y, s=24, color=C_LOY, zorder=3)
    else:
        axC.scatter(max(r.raw, 0) * 100, y, s=24, color="0.6", zorder=3)
    axC.text(108, y, f"{r.tf:.2f}", va="center", fontsize=6.5, color="0.4")
axC.text(108, -0.25, "TFx", fontsize=6.5, color="0.4", va="center", fontweight="bold")
axC.set_ylim(len(order), 0); axC.axhline(npos, color="0.15", lw=1, xmax=0.9)
axC.axvline(CUT * 100, color="0.5", lw=0.7, ls="--")
axC.set_xlim(-3, 116); axC.set_xticks([0, 25, 50, 75, 100]); axC.set_xlabel("Cells with chrY loss (%)")
axC.set_yticks(np.arange(len(order)) + 0.5); axC.set_yticklabels(["USC-" + l + ("*" if s in lowconf else "") for l, s in zip(lab, order)], fontsize=7); axC.tick_params(axis="y", length=0); axC.spines["left"].set_visible(False)
axC.text(-0.25, npos / 2, "LOY+", transform=axC.get_yaxis_transform(), rotation=90, ha="center", va="center", fontsize=7.5, fontweight="bold")
axC.text(-0.25, npos + (len(order) - npos) / 2, "LOY−", transform=axC.get_yaxis_transform(), rotation=90, ha="center", va="center", fontsize=7.5, fontweight="bold")
axC.legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec=C_LOY, ms=5, label="Whole specimen"),
                    Line2D([], [], marker="o", ls="", color=C_LOY, ms=5, label="Tumor cells (TFx-adjusted)"),
                    Line2D([], [], color=C_LOY, alpha=0.3, lw=5, label="Range across TFx scenarios"),
                    Line2D([], [], color="0.5", lw=0.7, ls="--", label="10% threshold")],
           frameon=False, fontsize=6.2, loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=2, handlelength=1.6, columnspacing=1.0)
axC.set_title("B", loc="left", fontsize=9, fontweight="bold", x=-0.3)
fig.savefig(FIGS / "main" / "Figure3.png", dpi=600, bbox_inches="tight"); fig.savefig(FIGS / "main" / "Figure3.pdf", bbox_inches="tight")
