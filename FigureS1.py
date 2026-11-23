import pandas as pd, numpy as np, matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from pathlib import Path
ROOT = Path("/Users/donny/LOY_UM")  # local clone of this repository; edit if cloned elsewhere
DATA, RESULTS, FIGS = ROOT / "data", ROOT / "results", ROOT / "figures"
for _d in (RESULTS / "Source_Data", FIGS / "main", FIGS / "supplementary"): _d.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.5})
X = DATA / "DataS1_segmented_copy_number_and_dPCR.xlsx"
M = RESULTS / "Source_Data" / "LOY_UM_MasterData.xlsx"
t = pd.read_excel(X, "tumor_TFx_adjusted_seg"); e = pd.read_excel(X, "excluded_tumor_seg_027_051"); b = pd.read_excel(X, "blood_seg")
keep = ~((t.chrom == 24) & ((t.chrompos < 2.8e6) | (t.chrompos > 26.6e6)))
t, b = t[keep].reset_index(drop=True), b[keep].reset_index(drop=True)
ts = pd.read_excel(M, "tumor_sWGS_summary").set_index("Patient"); bs = pd.read_excel(M, "blood_sWGS_summary").set_index("Patient")
order = ["140", "047", "045", "033", "064", "071", "068", "048", "005"]
chrom = t.chrom.values; n = len(chrom)
exp = np.where(chrom >= 23, 1.0, 2.0)
bd = np.r_[0, np.where(np.diff(chrom) != 0)[0] + 1, n]
iX = np.where(chrom == 23)[0]; iY = np.where(chrom == 24)[0]
cmap = LinearSegmentedColormap.from_list("cn", ["#1F3B73", "#6E8FC2", "#FFFFFF", "#F2B27A", "#C8641E"]); norm = TwoSlopeNorm(0, -2, 1)
RED = "#B03030"; C_A, C_X, C_Y = "#3E4C5A", "#5F84A8", "#C0392B"
fmt = lambda v: (f"{round(v):d}".replace("-", "−") if round(v) != 0 else "0") + "%"
lg = lambda v: np.log2(np.clip(2 * v / exp, 0.03, None))
T = np.array([lg(t[f"UM_{p}"].values) for p in order]); B = None
show = {1: "1", 3: "3", 6: "6", 8: "8", 12: "12", 16: "16", 20: "20", 23: "X", 24: "Y"}

def heat_pair(axG, axZ, Mx, title, labels):
    axG.imshow(Mx, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest", extent=(0, n, len(order), 0))
    for v in bd[1:-1]: axG.axvline(v, color="0.6", lw=0.25)
    for k in (7,): axG.axhline(k, color="0.1", lw=0.9)
    cs = list(pd.unique(chrom)); axG.set_xticks((bd[:-1] + bd[1:]) / 2); axG.set_xticklabels([show.get(c, "") for c in cs], fontsize=6.3)
    axG.set_yticks(np.arange(len(order)) + 0.5); axG.set_yticklabels([f"USC-{p}" + ("*" if p in ("033", "048") else "") for p in order], fontsize=6.5)
    axG.tick_params(length=0); [s.set_visible(False) for s in axG.spines.values()]
    axG.set_title(title, loc="left", fontsize=8, fontweight="bold")
    # zoom: X and Y with equal display width
    ix = iX[np.linspace(0, len(iX) - 1, 40).round().astype(int)]; iy = iY[np.linspace(0, len(iY) - 1, 40).round().astype(int)]
    Z = np.hstack([Mx[:, ix], Mx[:, iy]])
    axZ.imshow(Z, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest", extent=(0, 80, len(order), 0))
    axZ.axvline(40, color="0.3", lw=0.6)
    for k in (7,): axZ.axhline(k, color="0.1", lw=0.9)
    axZ.set_xticks([20, 60]); axZ.set_xticklabels(["X", "Y"], fontsize=6.8, fontweight="bold"); axZ.set_yticks([])
    axZ.tick_params(length=0)
    for s in axZ.spines.values(): s.set_edgecolor(RED); s.set_linewidth(0.9)
    for i, (txt, col) in enumerate(labels): axZ.text(86, i + 0.5, txt, va="center", fontsize=6.3, color=col)
    for (a, z, lab) in [(0, 7, "LOY+"), (7, 9, "LOY−")]:
        axG.text(-0.12, (a + z) / 2, lab, transform=axG.get_yaxis_transform(), rotation=90, ha="center", va="center", fontsize=6.8, fontweight="bold")


tfa = pd.read_excel(M, "Fig2C_LOY_TFx") if "Fig2C_LOY_TFx" in pd.ExcelFile(M).sheet_names else None
labs = [(fmt(ts.loc['USC-'+p, 'raw_LOY_fraction_pct']), RED if ts.loc['USC-'+p, 'raw_LOY_fraction_pct'] >= 10 else "0.3") for p in order]
fig = plt.figure(figsize=(7.2, 3.2))
heat_pair(fig.add_axes([0.13, 0.3, 0.67, 0.62]), fig.add_axes([0.81, 0.3, 0.045, 0.62]), T, "", labs)
fig.text(0.862, 0.935, "LOY", fontsize=6.5, fontweight="bold", color="0.3")
cax = fig.add_axes([0.32, 0.12, 0.3, 0.028]); sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
cb.set_ticks([-2, -1, 0, 1]); cb.ax.tick_params(labelsize=5.8, length=2); cb.outline.set_linewidth(0.4)
cb.set_label("log2 tumor-cell copy state relative to the normal male state (TFx-adjusted)", fontsize=6)
fig.savefig(FIGS / "supplementary" / "FigureS1.png", dpi=600); fig.savefig(FIGS / "supplementary" / "FigureS1.pdf")
