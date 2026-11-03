# Quarantine Games on Networks

Simulation, analysis and retained scientific data for *Quarantine Games on
Networks: Beyond Epidemic Establishment* by Marco Antonio Amaral.

**Version:** `rc9.5-p10.4-20261103` — reproducibility prerelease accompanying the RC9.5 P10.4
manuscript candidate. [Download the matching data](https://github.com/AmaralMA42/quarantine_games_network/releases/download/rc9.5-p10.4-20261103/quarantine_games_network_rc9.5_p10.4_data_20261103.zip)
(156.5 MB) or visit the [versioned release](https://github.com/AmaralMA42/quarantine_games_network/releases/tag/rc9.5-p10.4-20261103).

## Reproduce the figures

Use Python 3.12. Install the pinned dependencies, download the data ZIP without
extracting it, and run from the repository directory:

```text
python -m pip install -r requirements.txt
python main.py prepare-data --archive /path/to/quarantine_games_network_rc9.5_p10.4_data_20261103.zip --check-only
python main.py prepare-data --archive /path/to/quarantine_games_network_rc9.5_p10.4_data_20261103.zip
python main.py reproduce-published
```

On Windows, quote paths containing spaces. Preparation validates SHA-256 checksums
for the archive, all 1,366 manifest entries and 1,313 retained file bindings. It
preserves the checked-out source files and stops on conflicting scientific inputs.
The commands use stored results and do not rerun epidemic campaigns.

The reconstruction creates **17 PDF assets: 11 main and five supplementary
figures**, with two assets for Fig. 5. Main figures are written to
`reproduction/main/manuscript/figures/`; supplementary figures to
`reproduction/si/supplement/figures/`. Numerical checks and receipts are in `qa/`.
The [figure index](docs/FIGURE_INDEX.csv) locates inputs and active functions.

## One model implementation

`main.py` is the entry point. A single C++17 CTMC kernel in
`engine/src/fpm_network_sim.cpp` implements the stochastic rules for network and
finite well-mixed experiments. Python modules manage parameters, execution,
graph-first means, measured-cell matrices and scientific presentation. There is
no simulator selected by figure number. Archived renderers in the data package
record provenance and are not imported or executed by the active pipeline.

```text
python main.py --help
python main.py build-engine --compiler PATH_TO_CXX_COMPILER
python main.py validate --design EXPERIMENTS.tsv
python main.py simulate --design EXPERIMENTS.tsv --output NEW_RESULTS.tsv
```

A C++17 compiler is required for new experiments, not for reconstruction from
stored data. See [VARIABLES.md](VARIABLES.md) for the run specification and
[reproducibility notes](docs/REPRODUCIBILITY.md) for analysis commands and scope.

## Data and interpretation

The data archive includes measured results, graph edge lists, derived tables,
sampling metadata, reference figures and provenance. Missing heatmap cells are
left empty. Parameter guides join measured coordinates without fitted smoothing.
Network uncertainty retains the documented graph hierarchy and conditioning.
The archive is a licensed derivative of an immutable retained scientific package;
all original 1,313 bindings are unchanged. [Data record](DATA_BUNDLE.json).

Historical file prefixes and the Python namespace `fpm` are retained for
compatibility. This repository concerns this paper; the FPM framework is separate.

## License and citation

Original code and software documentation: [MIT](LICENSE), preserving the license
chosen when this repository was created. Original scientific data, tables and
figures: [CC BY 4.0](LICENSES/CC-BY-4.0.txt). See [license scope](LICENSING.md).
Please cite this version using [CITATION.cff](CITATION.cff). No DOI is claimed.

This prerelease preserves the manuscript's candidate status. It does not represent
journal acceptance or automatic promotion of the scientific version. The original
deterministic-model repository remains at
[AmaralMA42/quarantine_game](https://github.com/AmaralMA42/quarantine_game).
