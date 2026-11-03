# Reproducibility and analysis

## General commands

These commands consume user-supplied run tables and grids described in
`VARIABLES.md`. Output paths for new experiments must be unused.

```text
python main.py plan --template TEMPLATE.tsv --grid GRID.json --output NEW_DESIGN.tsv
python main.py validate --design NEW_DESIGN.tsv
python main.py simulate --design NEW_DESIGN.tsv --output NEW_RESULTS.tsv
python main.py averages --input NEW_RESULTS.tsv --output MEANS.csv --groups topology beta_n behavior_mode --values ever_infected_count
python main.py heatmap --input ONE_NETWORK_ONE_OBSERVABLE.csv --output MATRIX.csv --x beta_n --y delta --value mean
python main.py render --group all
python verify_inputs.py --include-generated
```

Means first average within observed graphs and then weight graphs equally. The
generic standard error is not a confidence interval. Matrices preserve missing
coordinates rather than interpolating. `render` draws the stored figures;
`reproduce-published` additionally runs the retained numerical and provenance checks.

## Scope of the retained evidence

The archive contains the stored inputs used to regenerate the selected figures;
it is not a claim that every historical campaign has been rerun or that every
simulator mode has received exhaustive boundary testing. A fresh reconstruction
was tested on Windows with Python 3.12 and the pinned dependency versions.
Cross-platform byte identity is not promised; compare measured values and plots
when using a different numerical, font or graphics environment.

Figure 3 A/B and Fig. 10 retain intervals conditional on observed graph samples.
Fig. 3 C/D and Fig. 11 have no error bars. Landmark selections retain their original
scope. The S5 pilot uses four graphs and remains exploratory. Recovered S4 records
are terminal outcomes, not 832 complete time series. Some historical graph-generator
seeds are unavailable and are not reconstructed by assumption.

The main simulations with initial infected fraction 0.01 and the supplementary
single-introduction experiments are distinct designs. Descriptive ordering within
degree families is not a universal ranking or a linear law. Independent comparison
samples and retained discrepancies remain separate from the published estimates.

## Files and checksums

`SOURCE_MANIFEST_SHA256.csv` covers every tracked source/document file except
itself. `DATA_BUNDLE.json` identifies the exact data asset. Its internal manifest
covers every archive member except itself; `provenance/SELECTED_INPUTS.csv`
additionally preserves 1,313 original bindings. The data attachment command checks
these before writing and does not replace repository sources.

Generated scientific inputs and outputs are ignored by Git. The data asset also
retains dated provenance notices; current public licensing is defined by
`LICENSING.md`. No editorial conversations, author handoff package, publisher
classes, compiled executables or FPM framework distribution are included.
