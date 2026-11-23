import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import binomtest, fisher_exact

from pathlib import Path
ROOT = Path("/Users/donny/LOY_UM")  # local clone of this repository; edit if cloned elsewhere
DATA, RESULTS, FIGS = ROOT / "data", ROOT / "results", ROOT / "figures"
for _d in (RESULTS / "Source_Data", FIGS / "main", FIGS / "supplementary"): _d.mkdir(parents=True, exist_ok=True)
def wilson(k, n, z=1.959964):
    p = k / n; den = 1 + z*z/n
    c = (p + z*z/(2*n)) / den
    h = z * ((p*(1-p)/n + z*z/(4*n*n)) ** 0.5) / den
    return c - h, c + h

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6})

h = pd.read_csv(DATA / "horsman1993_recode.csv"); h["patient"] = h.patient.astype(str); h["study"] = "Horsman 1993"
p = pd.read_csv(DATA / "prescher1990_recode.csv"); p["study"] = "Prescher 1990"
cols = ["study", "patient", "age", "sex", "mono3", "gain8q", "sex_chr_loss", "sex_chr_loss_type"]
d = pd.concat([h[cols], p[cols]], ignore_index=True)
assert p.mono3.sum() == 6 and p.gain8q.sum() == 8 and (p.sex_chr_loss_type == "LOY").sum() == 4 and (p.sex_chr_loss_type == "LOX").sum() == 1
d.to_csv(DATA / "UM_karyotype_pooled_Horsman1993_Prescher1990.csv", index=False)

C_LOY, C_LOX, C_M3, C_8Q, C_NO = "#B5545C", "#DE9A9F", "#4A5A6A", "#E08A3C", "#E6E6E6"

fig = plt.figure(figsize=(7.4, 5.8))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.72], width_ratios=[1.1, 1], hspace=0.62, wspace=0.3)

# A: per study + pooled
axA = fig.add_subplot(gs[0, 0])
sets = [("Horsman\n1993 (n=20)", d[d.study == "Horsman 1993"]), ("Prescher\n1990 (n=14)", d[d.study == "Prescher 1990"]), ("Pooled\n(n=34)", d)]
w = 0.34
for i, (lab, g) in enumerate(sets):
    n = len(g); k3 = g.mono3.sum(); ky = (g.sex_chr_loss_type == "LOY").sum(); kx = (g.sex_chr_loss_type == "LOX").sum()
    lo, hi = wilson(k3, n)
    axA.bar(i - w/2, k3/n*100, w, color=C_M3, yerr=[[k3/n*100 - lo*100], [hi*100 - k3/n*100]], ecolor="0.25", capsize=2, error_kw=dict(lw=0.6))
    axA.bar(i + w/2, ky/n*100, w, color=C_LOY)
    axA.bar(i + w/2, kx/n*100, w, bottom=ky/n*100, color=C_LOX)
    ks = ky + kx; lo, hi = wilson(ks, n)
    axA.errorbar(i + w/2, ks/n*100, yerr=[[ks/n*100 - lo*100], [hi*100 - ks/n*100]], color="0.25", capsize=2, lw=0.6)
    axA.text(i - w/2, 2, f"{k3}", ha="center", fontsize=6.5, color="white")
    axA.text(i + w/2, 2, f"{ks}", ha="center", fontsize=6.5, color="white")
axA.set_xticks(range(3)); axA.set_xticklabels([s[0] for s in sets])
axA.set_ylim(0, 100); axA.set_ylabel("% of abnormal tumors (95% CI)")
axA.legend(handles=[Patch(color=C_M3, label="−3"), Patch(color=C_LOY, label="−Y (♂)"), Patch(color=C_LOX, label="−X (♀)")],
           frameon=False, fontsize=7, ncol=3, loc="upper center", handlelength=1)
axA.set_title("A", loc="left", fontsize=9, fontweight="bold", x=-0.12)

# B: pooled, by sex, eligible denominators
axB = fig.add_subplot(gs[0, 1])
grp = [("Male\n(n=20)", d[d.sex == "M"], C_LOY), ("Female\n(n=14)", d[d.sex == "F"], C_LOX)]
for i, (lab, g, col) in enumerate(grp):
    n = len(g)
    for off, k, c in [(-w/2, g.mono3.sum(), C_M3), (w/2, g.sex_chr_loss.sum(), col)]:
        lo, hi = wilson(k, n)
        axB.bar(i + off, k/n*100, w, color=c, yerr=[[k/n*100 - lo*100], [hi*100 - k/n*100]], ecolor="0.25", capsize=2, error_kw=dict(lw=0.6))
        axB.text(i + off, 2, f"{k}/{n}", ha="center", fontsize=6.5, color="white")
    # per-study points
    for s, mk in [("Horsman 1993", "o"), ("Prescher 1990", "^")]:
        gs_ = g[g.study == s]
        for off, k in [(-w/2, gs_.mono3.sum()), (w/2, gs_.sex_chr_loss.sum())]:
            axB.scatter(i + off, k/len(gs_)*100, marker=mk, s=12, facecolor="white", edgecolor="0.15", lw=0.6, zorder=5)
axB.set_xticks(range(2)); axB.set_xticklabels([g[0] for g in grp])
axB.set_ylim(0, 100); axB.set_ylabel("% of tumors (95% CI)")
axB.legend(handles=[plt.Line2D([], [], marker="o", ls="", mfc="white", mec="0.15", ms=4, label="Horsman"),
                    plt.Line2D([], [], marker="^", ls="", mfc="white", mec="0.15", ms=4, label="Prescher")],
           frameon=False, fontsize=7, ncol=2, loc="upper center", handlelength=1)
axB.set_title("B", loc="left", fontsize=9, fontweight="bold", x=-0.12)

# C: oncoprint
axC = fig.add_subplot(gs[1, :])
parts = []
for s in ["Horsman 1993", "Prescher 1990"]:
    for sx in ["M", "F"]:
        parts.append(d[(d.study == s) & (d.sex == sx)].sort_values(["sex_chr_loss", "mono3", "gain8q"], ascending=False))
order = pd.concat(parts)
rows = [("Sex chr loss", "sex_chr_loss"), ("−3", "mono3"), ("+8/8q", "gain8q")]
N = len(order)
for j, (_, r) in enumerate(order.iterrows()):
    for i, (lab, col) in enumerate(rows):
        if col == "sex_chr_loss":
            c = (C_LOY if r.sex_chr_loss_type == "LOY" else C_LOX) if r[col] else C_NO
        else:
            c = (C_M3 if col == "mono3" else C_8Q) if r[col] else C_NO
        axC.add_patch(plt.Rectangle((j + 0.07, -i - 0.88), 0.86, 0.76, color=c, lw=0))
    axC.add_patch(plt.Rectangle((j + 0.07, 0.25), 0.86, 0.45, color="#5B7FA6" if r.sex == "M" else "#C48BB0", lw=0))
axC.axvline(20, color="0.3", lw=0.8)
axC.set_xlim(0, N); axC.set_ylim(-len(rows) - 0.05, 0.8)
axC.set_yticks([0.47] + [-i - 0.5 for i in range(len(rows))]); axC.set_yticklabels(["Sex"] + [r[0] for r in rows])
axC.set_xticks(np.arange(N) + 0.5)
axC.set_xticklabels([str(x).replace("AM", "") for x in order.patient], fontsize=5.5, rotation=90)
axC.set_xlabel("Horsman 1993 (patient no.)                                                        Prescher 1990 (AM no.)", fontsize=7)
for s_ in ["left", "bottom"]: axC.spines[s_].set_visible(False)
axC.tick_params(length=0)
axC.legend(handles=[Patch(color="#5B7FA6", label="Male"), Patch(color="#C48BB0", label="Female"),
                    Patch(color=C_LOY, label="−Y"), Patch(color=C_LOX, label="−X"), Patch(color=C_NO, label="Absent")],
           frameon=False, fontsize=7, ncol=5, loc="lower center", bbox_to_anchor=(0.5, 1.0), handlelength=1)
axC.set_title("C", loc="left", fontsize=9, fontweight="bold", x=-0.12, pad=16)

fig.savefig(FIGS / "supplementary" / "FigureS4.png", dpi=600, bbox_inches="tight")
fig.savefig(FIGS / "supplementary" / "FigureS4.pdf", bbox_inches="tight")

for name, g in [("Prescher", d[d.study == "Prescher 1990"]), ("Pooled", d)]:
    both = ((g.mono3 == 1) & (g.sex_chr_loss == 1)).sum(); m3 = ((g.mono3 == 1) & (g.sex_chr_loss == 0)).sum()
    sc = ((g.mono3 == 0) & (g.sex_chr_loss == 1)).sum(); ne = ((g.mono3 == 0) & (g.sex_chr_loss == 0)).sum()
    print(name, "both/m3only/sconly/neither", both, m3, sc, ne, "McNemar p", round(binomtest(sc, m3 + sc).pvalue, 3),
          "Fisher", fisher_exact([[both, m3], [sc, ne]]))
    for sx in ["M", "F"]:
        gg = g[g.sex == sx]; print("  ", sx, "-3", gg.mono3.sum(), "/", len(gg), "SCL", gg.sex_chr_loss.sum())
pm = p[p.sex == "M"]
print("Prescher male ages LOY+", sorted(pm[pm.sex_chr_loss == 1].age), "LOY-", sorted(pm[pm.sex_chr_loss == 0].age))
