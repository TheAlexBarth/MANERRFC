# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Analysis code for a manuscript on microplankton trophic roles across wet/dry regime cycles, using Mission
Aransas NERR / SWMP monitoring data. All processing is in R; Bayesian models are written in Stan and fit
via `cmdstanr`. This is an analysis pipeline, not an R package — there is no `DESCRIPTION`, `renv.lock`,
test suite, or lint config. `! Warning` in README.md: analyses are ongoing/incomplete.

## Running scripts

There is no build system — run scripts directly, e.g. `Rscript R/fig-02-environ_ts.R` or source them in
RStudio (the repo has an `.Rproj`). Every script starts with `rm(list = ls())` and expects to be run
standalone in a fresh session with the working directory at the repo root (paths are written as `./R/...`,
`./data/...`, `./stan/...`). Scripts are not designed to be `source()`-chained into one another except via
the shared `R/utils*.R` files.

Packages are installed ad hoc into the local R library (no renv/lockfile). `EcotaxaTools` is a
non-CRAN/custom package this project depends on for Ecotaxa-format data handling and posterior-array
utilities (`R/pipe-00-import_etx.R`, `R/utils-posterior_specific.R`) — it must already be installed
locally.

## Script naming convention (see README.md)

- `pipe-##[letter]-*.R` — import and reformat one external raw data source, save an RDS/RDS-like file into
  `data/` named with the same numeric prefix (e.g. `pipe-01a-import_swmp_environ.R` → `data/01a-swmp_wq_data.rds`).
- `##-*.R` (no prefix word, e.g. `03-gen_taxa_pred-mod.R`, `04-troph_regression.R`) — analysis/model-fitting
  scripts. These load merged data, build `cmdstanr` models from `stan/*.stan`, and fit them.
- `fig-##-*.R` — generate a manuscript figure, written to `output/` as a PDF.
- `s-*.R` — supplementary/summary scripts and plots, also written to `output/`.
- `xx-*.csv` in `data/` — static external reference data that's small enough to live in the repo directly
  (site metadata, taxonomy maps, drought records) rather than being produced by a `pipe-` script.
- `utils*.R` — shared helpers sourced by other scripts (`source('./R/utils.R')`), never run directly.
- Scripts use `# region \-` / `# endregion` and `# MARK:` comments as foldable section markers (RStudio
  code-folding convention) — preserve this style when adding sections to existing files.

## Data pipeline order

1. `pipe-00-import_etx.R` — reads raw Ecotaxa plankton-imaging data and `data/xx-taxo_map.csv` →
   `data/00-ecotaxa_full.rds`.
2. `pipe-01a` (SWMP water quality), `pipe-01b` (size-fractionated chlorophyll), `pipe-01c` (wind, scraped
   from NOAA NDBC), `pipe-01d` (river regime), `pipe-01e` (silicate, SiOH) — each independent, →
   `data/01[a-e]-*.rds`.
3. `pipe-02-data_merger.R` — joins `00` + `01a/b/c/e` outputs into `data/02-full_merged.rds` (a list with
   `conc` and `indv` tibbles), which is the primary input for downstream analysis/model scripts.
4. `03-*` / `04-*` — fit Stan models (via `R/utils-mod_making.R::fit_mod()` and friends) → save posterior
   output to `data/03-post_chl.RDS`, `data/04-post_micro.RDS`, etc.
5. `fig-*` / `s-*` — consume the merged and posterior data to produce figures/tables in `output/`.

**Important:** `pipe-00` and `pipe-01a`/`pipe-01e` (and `fig-01-map.R`'s basemap layer) read raw source
files from a hardcoded local Box Cloud Storage path
(`~/Library/CloudStorage/Box-Box/TGCRC Plankton Food Webs/Data/...`) that only exists on the author's
machine — these scripts explicitly note "will need local adjustment" and cannot be rerun without that
access. The `data/*.rds` outputs of the `pipe-` stage are already committed to the repo, so
`03-`/`04-`/`fig-`/`s-` scripts can be run/modified without needing the raw data or Box access.

## Nutrient covariates vs. nutrient ratios

Regression covariates (`03-size_frac_chla_reg.R`, `04-troph_regression.R`) always use raw nutrient
*concentrations* (P, NH4, N, SiOH — log-transformed, then scaled), never molar ratios. Molar N:P / N:Si
stoichiometry (via `compute_molar_ratios()` in `R/utils.R`) is a *reporting-only* quantity: it's visualized
and trend-tested in `fig-02-environ_ts.R` (panels I/J, plus `output/fig02-nutrient_ratio_trends.csv` for
the per-site linear trend stats), not used as a model predictor. If you're extending the nutrient
regressions, keep this split — don't substitute a ratio term into `preds`/`select()` in the `03-`/`04-`
scripts.

## Modeling code structure

- `R/utils-mod_making.R` — data-prep helpers (`prep_component_data()`, `prep_component_data_sites()`) that
  build model-ready tibbles from the `etx` merged-data object, plus `fit_mod()` which assembles the
  `cmdstanr` data list and calls `$sample()`.
- `R/utils-posterior_specific.R` — posterior-processing/plotting helpers built on top of fitted draws:
  `post_predict()` (evaluate a formula over posterior draws + new predictor data), `summarize_pred()`
  (draws matrix → mean/median/HDI or quantile interval summary), `post_arr_to_df()`, `make_marg_data()` /
  `make_marg_plot()` (marginal-effect simulation and plotting), `make_ts_data()` / `make_ts_plot()`
  (time-series prediction and plotting).
- `stan/*.stan` — model definitions, compiled on demand via `cmdstanr::cmdstan_model()`. The sibling
  extensionless directories (`stan/bmass_mod`, `stan/count_mod`, etc.) are cmdstan-compiled binaries, not
  source — don't hand-edit them.
- Shared plotting constants (site/taxon/size factor levels and colorblind-friendly palettes: `site_cols`,
  `troph_cols`, `size_cols`, `regime_cols`, `site_factors`, `troph_factors`, `size_factors`) live in
  `R/utils.R` and are reused across all `fig-`/`s-` scripts for consistent figure styling.
- Figure PDFs are consistently saved at 85×85mm, 600dpi (`ggsave(..., width = 85, height = 85, units = 'mm', dpi = 600)`) — match this when adding new figure scripts unless the figure needs a different aspect ratio.
