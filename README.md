# LOY_UM

Analysis and plotting code for:

> Liang D, Singh M, Feng Z, Romero I, Galarneau H, Peng C-C, Kuhn P, Berry JL, Xu L. **Tumor and host loss of chromosome Y can be discordant in uveal melanoma.** *Scientific Reports* (submitted).

The code quantifies loss of chromosome Y (LOY) in intraocular specimens (aqueous humor and tumor wash) and matched blood from men with uveal melanoma using shallow whole-genome sequencing (sWGS), and measures mosaic LOY (mLOY) in blood by digital PCR (dPCR). It reproduces the figures, supplementary figures and source data tables of the article.

## What the code does

| Step | Where |
|---|---|
| Chromosome Y dosage relative to chromosome X (chrY 2.8–26.6 Mb vs. all chrX bins) | `build_master_data.py` |
| Empirical no-LOY baseline from matched blood without LOY | `build_master_data.py` |
| Raw LOY fraction, `L = 1 − (chrY:chrX) ÷ baseline`, and calls at the 10% threshold | `build_master_data.py` |
| One-sided normal tolerance factors (noncentral *t*) for the threshold calibration | `build_master_data.py`, `Figure4.py` |
| TFx-adjusted LOY fraction, `f = L ÷ TFx`, under primary and alternative TFx scenarios | `build_master_data.py` |
| Arm-level copy-number calls relative to tumor fraction | `build_master_data.py`, `Figure3.py` |
| dPCR mLOY, `1 − AMELY ÷ AMELX`, delta-method 95% CIs and assay calling threshold | `build_master_data.py` |
| Statistics (Wilcoxon signed-rank, Spearman, Clopper–Pearson, Bland–Altman) | `build_master_data.py`, `update_master_v5.py`, figure scripts |

Upstream processing (alignment with BWA-MEM, binning into 5000 variable-length bins, GC correction and segmentation with DNAcopy) follows Baslan et al. (*Nat. Protoc.* 2012) and is not included here; the scripts start from segmented copy ratios. Tumor fraction values from the clonal anchor event are read from the `Supp_TF_sensitivity` sheet of `seg.xlsx`.

## Scripts and manuscript figures

| Script | Manuscript item |
|---|---|
| `Figure2.py` | Fig. 2: genome-wide copy number of tumors and matched blood |
| `Figure4.py` | Fig. 3: tumor vs. matched blood LOY, threshold calibration, threshold sweep |
| `Figure5.py` | Fig. 4: host mLOY by dPCR, age, cross-platform concordance |
| `Figure3.py` | Fig. 5: LOY and arm-level alteration frequency, TFx-adjusted LOY |
| `FigureS2.py` | Supplementary Fig. S1: matched blood LOY under alternative baselines |
| `FigureS3_TableS3.py` | Supplementary Fig. S2 and Supplementary Table S1: sWGS vs. dPCR agreement |
| `FigureS1.py` | Supplementary Fig. S3: genome-wide TFx-adjusted tumor-cell copy state |
| `FigureS4.py` | Supplementary Fig. S4: sex chromosome loss in historical cytogenetic series |
| `build_master_data.py` | Source data workbook (`LOY_UM_MasterData.xlsx`), one sheet per panel |
| `update_master_v5.py` | Adds the cohort table, figure index and recomputed key numbers to the workbook |

Script names follow an earlier figure order; the table above gives the corresponding items in the article. Figure 1 is a schematic.

## Folder layout

All paths are set relative to one project folder, defined at the top of each script:

```python
ROOT = Path("/Users/donny/LOY_UM")
```

Change this line if the repository is cloned elsewhere. The scripts expect:

```
LOY_UM/
├── data/                      input tables (not tracked; see below)
│   ├── seg.xlsx                         tumor segmented copy ratios (raw, TFx-adjusted) and TFx scenarios
│   ├── blood.csv                        matched blood segmented copy ratios (11 patients)
│   ├── seg_027_051.csv                  segments of the two ploidy-indeterminate tumors
│   ├── DataS1_segmented_copy_number_and_dPCR.xlsx
│   ├── LOY_PCR.csv                      QIAcuity dPCR export
│   ├── dpcr_clean.csv                   per-sample dPCR summary
│   ├── FigS1_summary.csv                chrY and chrX dosage of the ploidy-indeterminate tumors
│   ├── horsman1993_recode.csv           recoded karyotypes, Horsman & White 1993
│   └── prescher1990_recode.csv          recoded karyotypes, Prescher et al. 1990
├── results/                   source data workbooks (created by the scripts)
└── figures/
    ├── main/                  (created by the scripts)
    └── supplementary/         (created by the scripts)
```

`data/`, `results/` and `figures/` are excluded from version control by `.gitignore`.

## Data

Segmented copy-number data and digital PCR measurements are provided in Supplementary Data S1 of the article. Raw sequencing data contain potentially identifiable genomic information and are available under controlled access as described in the article's Data availability statement. The historical karyotype tables were recoded from the published case listings of Prescher et al. (*J. Natl. Cancer Inst.* 1990) and Horsman & White (*Cancer* 1993).

## Requirements

Python 3 with:

```
pip install pandas numpy scipy matplotlib openpyxl
```

## How to run

From the repository folder, in this order:

```bash
python FigureS4.py            # also writes the pooled historical karyotype table used below
python build_master_data.py   # builds results/Source_Data/LOY_UM_MasterData.xlsx
python update_master_v5.py    # adds cohort, index and key-number sheets
python Figure2.py
python Figure3.py
python Figure4.py
python Figure5.py
python FigureS1.py
python FigureS2.py
python FigureS3_TableS3.py
```

Figures are written as PNG (600 dpi) and PDF to `figures/main/` and `figures/supplementary/`.

## Contact

Donny Liang (dliang@chla.usc.edu) or Liya Xu (lixu@chla.usc.edu), Children's Hospital Los Angeles and University of Southern California.
