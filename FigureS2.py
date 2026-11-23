import pandas as pd, numpy as np, matplotlib.pyplot as plt
from pathlib import Path
ROOT = Path("/Users/donny/LOY_UM")  # local clone of this repository; edit if cloned elsewhere
DATA, RESULTS, FIGS = ROOT / "data", ROOT / "results", ROOT / "figures"
for _d in (RESULTS / "Source_Data", FIGS / "main", FIGS / "supplementary"): _d.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6})
d = pd.read_excel(RESULTS / "Source_Data" / "LOY_UM_MasterData.xlsx", "FigS2A_blood_baselines")
d = d.sort_values("LOY_pct_baseline_primary", ascending=False).reset_index(drop=True)
C = "#6F8BA4"
fig, ax = plt.subplots(figsize=(7.2, 2.8))
x = np.arange(len(d))
ax.scatter(x - 0.18, d.LOY_pct_baseline_primary, s=26, color=C, label="Primary baseline 0.952 (n = 10)", zorder=3)
ax.scatter(x, d["LOY_pct_baseline_0.936"], s=26, marker="s", facecolor="white", edgecolor=C, lw=1, label="Alternative baseline 0.936 (n = 3)", zorder=3)
ax.scatter(x + 0.18, d["LOY_pct_baseline_1.0"], s=28, marker="^", facecolor="white", edgecolor=C, lw=1, label="Theoretical baseline 1.0", zorder=3)
ax.axhline(10, color="0.45", lw=0.8, ls="--"); ax.axhline(0, color="0.85", lw=0.6)
ax.text(len(d) - 0.5, 10.4, "10% threshold", ha="right", va="bottom", fontsize=6.5, color="0.4")
ax.set_xticks(x); ax.set_xticklabels(["USC-" + str(p).replace("USC-", "") for p in d.Patient], fontsize=7)
ax.set_ylim(-8, 20); ax.set_ylabel("Blood LOY fraction (%)"); ax.set_xlim(-0.6, len(d) - 0.4)
ax.legend(frameon=False, fontsize=6.8, loc="upper right", ncol=1)
fig.savefig(FIGS / "supplementary" / "FigureS2.png", dpi=600, bbox_inches="tight"); fig.savefig(FIGS / "supplementary" / "FigureS2.pdf", bbox_inches="tight")
