# Quarantine Games on Networks

Simulation, analysis and retained data for *Quarantine Games on Networks:
Beyond Epidemic Establishment* by Marco Antonio Amaral.

**Version:** `rc9.5-p10.4-authorial-20261106`, accompanying the current authorial revision.
This repository concerns the paper; the FPM framework is separate. The manuscript
and journal-submission status are independent of the software release.

## Reproduce the current figures

Use Python 3.12 and the pinned dependencies. Download both archives without
extracting them: the unchanged [base data](https://github.com/AmaralMA42/quarantine_games_network/releases/download/rc9.5-p10.4-20261103/quarantine_games_network_rc9.5_p10.4_data_20261103.zip)
(156.5 MB) and the [current figure update](https://github.com/AmaralMA42/quarantine_games_network/releases/download/rc9.5-p10.4-authorial-20261106/quarantine_games_network_authorial_20261106_update.zip)
(about 1.55 MB). Then run from the repository directory:

```text
python -m pip install -r requirements.txt
python main.py prepare-data --archive PATH_TO_BASE.zip --update-archive PATH_TO_UPDATE.zip --check-only
python main.py prepare-data --archive PATH_TO_BASE.zip --update-archive PATH_TO_UPDATE.zip
python main.py reproduce-published
```

Quote paths containing spaces. Both archives are verified before attachment.
The base retains 1,366 manifest entries and 1,313 original bindings. The update
adds 51 current input/analysis files, 18 reference graphics and the current
figure index in a separate namespace; it does not overwrite repository sources
or historical scientific inputs. [Base record](DATA_BUNDLE.json),
[update record](ARTICLE_UPDATE.json).

The command generates **18 PDF assets for 11 main and six supplementary figures**;
Fig.5 uses two assets. Main outputs are in
`reproduction/main/manuscript/figures/`; SI outputs in
`reproduction/si/supplement/figures/`. QA receipts are in `qa/`.
The [figure index](docs/FIGURE_INDEX.csv) identifies the inputs and functions.
Current S6 contains relative-transmission bars. The earlier reach/protection
view and the historical four-panel/full-landmark layouts remain archived data,
outside the current selection.

## One model implementation

`main.py` is the entry point. The single C++17 CTMC kernel in
`engine/src/fpm_network_sim.cpp` implements the network and finite well-mixed
rules. Python modules manage parameters, execution, graph-first means,
measured-cell matrices and presentation. Shared visual modules compose all
current views; no simulator is selected by figure number. Archived renderer
scripts record provenance and are not imported or executed by the active pipeline.

```text
python main.py --help
python main.py build-engine --compiler PATH_TO_CXX_COMPILER
python main.py validate --design EXPERIMENTS.tsv
python main.py simulate --design EXPERIMENTS.tsv --output NEW_RESULTS.tsv
```

A compiler is required for new experiments, not for reconstruction from stored
data. See [VARIABLES.md](VARIABLES.md), [reproducibility notes](docs/REPRODUCIBILITY.md)
and the [scientific analysis manual](docs/SCIENTIFIC_METHODS.md).

## Evidence and interpretation

Missing heatmap cells remain empty; parameter guides join measured coordinates
without fitted smoothing. Means and retained intervals preserve their documented
graph hierarchy and conditioning. Fig.3 C/D and Fig.11 have no error bars.
The original observations and their provenance remain intact. Current generated
graphics are checked against their own index, separately from the historical
generated-output records; the original provenance ledger is not rewritten.

A clean reconstruction was verified on Windows with the recorded Python and
dependency versions. It does not rerun epidemic campaigns or validate every
historical simulation mode. Cross-platform byte identity is not promised.

## License and citation

Original code/software documentation: [MIT](LICENSE). Original scientific data,
tables and figures: [CC BY4.0](LICENSES/CC-BY-4.0.txt), in the scopes documented in
[LICENSING.md](LICENSING.md). Please cite this version through [CITATION.cff](CITATION.cff).
No article or data DOI is claimed. The [previous release](https://github.com/AmaralMA42/quarantine_games_network/releases/tag/rc9.5-p10.4-20261103)
remains immutable and supplies the base archive.

The original deterministic model is at
[AmaralMA42/quarantine_game](https://github.com/AmaralMA42/quarantine_game).
