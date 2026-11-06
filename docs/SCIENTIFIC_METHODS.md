# Scientific analysis protocols

Documentation supplement, 6 November 2026. Current authorial revision of
RC9.5 P10.4, with unchanged measured inputs and scientific procedures.

This manual retains analysis details condensed in the manuscript and SI.
The public source version is `rc9.5-p10.4-authorial-20261106`. The historical base
archive remains unchanged; `ARTICLE_UPDATE.json` identifies the small current
input/reference addition. No simulation rule, observation or interval endpoint
was changed for this release.

The entry point remains `main.py`, with shared dynamics and general analysis
modules. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for commands,
[VARIABLES.md](../VARIABLES.md) for the parameter contract and
[FIGURE_INDEX.csv](FIGURE_INDEX.csv) for the stored inputs and rendering functions.
Paths below refer to the attached data bundle unless explicitly described as code.
The manuscript and SI define the model; this manual provides the analysis recipe.

## 1. Dynamics and notation

The states are susceptible-normal, susceptible-protected, infected-normal,
infected-protected and recovered. Infection preserves strategy. Recovery is
absorbing and removes the individual from strategy updating. Infection, recovery
and imitation compete through their instantaneous rates in one continuous-time
Markov process. Gillespie waiting times are exponential with rate equal to the
sum of event propensities; event selection is proportional to the propensity.
Every scientific realization used here ends at infection-free absorption.

The standard settings are Omega = tau = gamma = 1, K = 0.1, beta_Q = 1,
alpha = 1/2, q_0 = 0.5 and i_0 = 0.01. The current mixed-strategy rate is
beta_A = (beta_N + beta_Q)/2. This differs from the extra factor a = 0.1 in the
mixed rate of the original deterministic publication. Infection and imitation
rates on networks are normalized by degree; finite well-mixed infection
propensities use N_0 - 1 possible partners. These normalizations are part of the
model, rather than adjustments to plotted data.

Global information uses population prevalence. Aligned local information uses
infected vertices at graph distances one through four, excluding the focal
vertex. A fixed permutation reassigns these local values in the shuffled control.
Each control calculates risk from its own evolving epidemic. Frozen controls
disable imitation and retain infection and recovery. Matched comparisons share
the graph and initial allocations; the retained paired designs also specify their
shared random seeds. A fixed strategy does not fix the active protected fraction,
because protected individuals can recover.

## 2. Observables and sampling

Reach is the fraction ever infected, including initial infections. Establishment
uses reach at least 0.05. Peak prevalence is the maximum infected fraction;
duration is time to absorption. Terminal active protection excludes recovered
individuals. Time-averaged protection is computed as the protection time integral
divided by that realization's duration, followed by averaging those ratios.
Cumulative net strategy balance is (N-to-Q events minus Q-to-N events)/N.
Means and recurrence frequencies include failed establishments.

Spatial results condition on a fixed periodic graph. For RRGs, average epidemics
within each graph first and weight the graph means equally. The sample size n
is the number of epidemic realizations per parameter condition, not the number
of independent graph realizations G. Different graph sets, population sizes,
initializations and detector protocols remain separate designs.

| Main figure or analysis | Epidemic realizations | RRG graphs / design |
|---|---|---|
| 3 A/B, initial and additional coordinates | 128 / 64 | G = 8; 16 / 8 per graph; graph sets distinct |
| 3 C/D; networks in 4 and 11 | 128 | G = 8; 16 per graph; N approximately 2^15 |
| 4, finite well-mixed | 16 | N = 4096 |
| 5 and 6 | 128 or 256 | G = 16; 8 or 16 per graph; N approximately 2^14 |
| 7, smaller and largest sizes | 64 / 256 | G = 4 / 8; 16 / 32 per graph; largest near 2^18 |
| 8, later-response landmarks | Initial Fig. 3 grid | Establishment reference from Fig. 7 A |
| 9, matched controls | 32 matched realizations | G = 8; four pairs per graph |
| 10 A / B | 64 / 32 | G = 4; 16 / 8 per graph |

Near nominal 2^15, actual populations are 33120 for 4.8.8, 33124 for Honeycomb,
32760 for Kagome, 32761 for Square and Triangular, and 32768 for RRGs. Fractions
use these actual counts. Main Fig. 4 A displays beta_N through 20;
main Fig. 4 B and Fig. 11 display measured values through 15; the stored
response and matched-degree series extend through 40. Fig. 3 C/D uses the initial
grid, not a pooled initial-plus-additional ensemble.

The heatmaps use a 28 by 20 measured grid per topology over beta_N from 1.5 to 15
and delta from 1 to 20. Additional realizations in selected cells produce the two
sample sizes. Retained intervals do not correct for that selection. Missing
coordinates are not interpolated, and scales are specific to each observable.

Graphs are reused across transmission values. The two degree-three contrasts in
Fig. 11 share the RRG3 reference, and the two degree-four contrasts share RRG4.
These are dependent comparisons, not five independent RRG experiments.

## 3. Probability intervals and mean contrasts

All probability intervals discussed here are pointwise, rather than simultaneous
over the grid, thresholds or networks. In Figs. 3 A/B and 10, the target conditions
on the observed graph or graph bank; no interval for a newly sampled graph
population is claimed.

For a fixed spatial graph or the finite well-mixed process, observed success
counts x of n use two-sided 95% Clopper-Pearson binomial intervals. The lower
endpoint is the 0.025 quantile of Beta(x, n-x+1), or zero when x = 0; the upper
endpoint is the 0.975 quantile of Beta(x+1, n-x), or one when x = n.

For the retained RRG probability comparisons, equal numbers of independent
epidemics per graph give the equal-graph mean q. Conservative heterogeneous
Bernoulli Chernoff-KL endpoints solve

    n * [q log(q/p) + (1-q) log((1-q)/(1-p))] = log(40).

There is one root on each side of q. For q = 0 the bounds are
[0, 1-exp(-log(40)/n)]; for q = 1 they are [exp(-log(40)/n), 1]. This construction
allows different probabilities on the fixed observed graphs. It assumes
conditional independence of the epidemic realizations. The implementation for
Fig. 3 is in `fpm/visuals/inference.py`; Fig. 10 consumes its retained conditional
interval table identified by the figure index.

Main Fig. 4 retains archived network-bootstrap endpoints
and finite well-mixed normal-approximation endpoints. Its input values are in
`reproduction/main/data/fig4_values.csv`, with source material identified by the
figure index. This editorial transfer supplies no new resampling algorithm or
seed for an incompletely documented historical calculation. Stored endpoint
reproduction and recovery of every historical bootstrap detail are distinct.
No missing recipe is filled in by inference from summary endpoints.

For main Fig. 9, form paired outcome differences at each fixed reference
condition. The retained procedure uses 1000 bootstrap replicates: resample
matched pairs on the spatial fixed graph; for RRGs, resample graphs first and
then matched pairs within the selected graphs. Intervals are pointwise 95% and
conditional on the selected reference coordinates. Supplementary Fig. S3 shows
the larger set of reference conditions under the same paired design. Positive
Frozen-minus-Adaptive, Local-minus-Global or Aligned-minus-Shuffled means larger
reach in the first condition. Opposite conventions must also reverse endpoints:
[lo, hi] becomes [-hi, -lo].

Fig. 11 shows the five spatial-minus-RRG mean contrasts without error bars, by
author choice. The full series and uncertainty estimates remain stored. Observed
sign changes in these means are not silently converted into resolved reversals.
Selected larger-population comparisons are separate experiments, reported in SI
section 5.

For independent validation of means, new experiments need not reuse historical
seeds. Comparisons should retain the estimand, initialization, graph hierarchy
and relevant conditions, with uncertainty appropriate to the newly sampled
graphs. Within a paired experiment, the pairing specified by that design remains
part of the analysis.

## 4. Peak detection and illustrative histories

For the selected recurrence-depth comparisons, the robust count is the
sixth-smallest of the peak counts from 27 configurations, capped at eight. Thus
at least 22 configurations must report m or more peaks for level m to be robust.
The configuration grid combines moving-average windows of 5, 9 and 17 samples;
relative prominence fractions 0.05, 0.10 and 0.20; and minimum separations of
12, 24 and 40 samples. With recording step 0.25, the separations are 3, 6 and
10 time units. Prominence is relative to each smoothed maximum and at least one
infected individual. Peaks must exceed max(1, N I-star), where
I-star = Omega/(delta beta_N) is the payoff-equality prevalence.

The 27 definitions are alternative measurements of a recorded history, not 27
independent epidemic realizations. P(n_p >= m) is estimated from the fraction of
epidemics meeting the criterion. Threshold counts for m = 2, 3, 4 and 5 use the
same epidemics and are dependent. The two regimes of main Fig. 10 are
(beta_N, delta) = (4, 25), n = 64, and (8, 10), n = 32. Both parameters change.

The dense recurrence maps retain their recorded detector settings and temporal
resolution. They are not retroactively treated as if they used an identical
selected-point measurement protocol. The mean peak-count map and the
probability of at least a second peak measure different aspects of recurrence.
Deep-trough return is a separate recorded observable. After the first reference
peak, prevalence must fall below 0.005 and later reach at least 0.01. The retained
summary field `reemergence_given_deep` conditions on histories with such a trough.
At (4, 25), this conditional return probability is 0.96875 for both RRG3 and RRG4,
whereas P(n_p >= 2) is 0.6875 and 0.921875 and P(n_p >= 3) is 0.046875 and 0.84375.
The return statistic must not be replaced by the robust second-peak frequency.
Its historical source is `M32EXT3_B4D25_MATCHED_MECHANISM_SUMMARY.csv` and the
threshold rule is in `m32ext3_mechanisms.py`, both retained in the public data
archive. This editorial clarification does not rerun the mechanism analysis.

The Fig. 2 examples were selected from six Square histories and twelve matched
RRG4 pairs on three graphs; all reached absorption. The displayed Square history
uses (4, 25). The displayed RRG4 adaptive/frozen pair uses (5, 16), with shared
graph, initial allocations and dynamics seed. The adaptive history has at least
four peaks under all 27 detector settings and the frozen history has one. A
separate frozen ensemble has no robust multipulse events in 288 realizations
over the studied conditions. Representative histories illustrate the process;
the ensemble supplies frequencies. The histories are actual simulations, not
curves drawn over an old image.

## 5. Crossovers, landmarks and deterministic endpoints

Establishment crossing coordinates beta_05, beta_50 and beta_95 correspond to
probabilities 0.05, 0.50 and 0.95 at the 5% reach threshold. Fig. 7 uses equal-weight
isotonic regression at smaller sizes and the cumulative maximum near N = 2^18,
followed by linear interpolation of crossings on these monotone curves.
The crossover width is beta_95 minus beta_05. Retained bootstrap bounds span
the central 95% of resampled widths for which both crossings are observed.

At the smallest population, valid fractions are 130/200 for RRG3 and 126/200 for
Triangular. Dotted bars and daggers mark these fractions below 95%; every other
displayed valid fraction is at least 95%. The intervals condition on successful
crossing observation and are not unconditional confidence bounds for all possible
crossovers. Stored width estimates and fractions are in
`reproduction/main/data/fig7_archived_widths.csv` and windows in `fig7_windows.csv`.

Main Fig. 8 displays three landmarks: the establishment midpoint near 2^18,
the post-establishment reach minimum and the time-averaged-protection maximum near
2^15. Both later markers select measured coordinates in 8 <= beta_N <= 40;
reach uses argmin of the mean and protection uses argmax of the mean of
per-realization time averages. A tie selects the smallest beta_N. These are
exactly the retained `selection_rule` fields, not a new optimization. The complete
48-row input table (six landmarks on eight networks), original intervals and
selection fields remain unchanged in `inputs/fig8_landmarks.csv`. The three
additional archived landmarks are the terminal-protection maximum, net-strategy-
balance maximum and rebound. The archived rebound is the first measured beta_N
whose lower reach bound exceeds
the upper bound at the selected reach minimum, using pointwise comparisons on
the observed graphs. The reach minimum for RRG6 is sensitive to the left edge of
its selection window. Some measured positions coincide. The 4.8.8/Honeycomb time-averaged-protection
maxima have distinct point estimates with overlapping location intervals. These are selected measured
coordinates, not a fitted universal degree law. Fig. S3's first subsequent
increase in mean reach is a different definition from this interval-separation
rebound. Fig. 8 retains its broken linear axis, 1 to 2.5 and 8 to 20, without bars.

For ODE duration in main Fig. 4 and the endpoint comparisons in Fig. S1,
evaluate the first time at which
infected density falls below 1/16384. Reach, duration and active protection use
that stopping time. If a local minimum drops below the threshold, first passage
switches to an earlier decline and the reported endpoint can jump despite the
continuous ODE solution. Fig. S1 changes the threshold to c/N for c = 0.5, 1 and 2;
the branch jump moves. These ODE endpoints and population sizes differ from the
finite well-mixed and network comparisons and must not be relabeled as identical
experiments.


The supplementary relative-landmark view (Fig. S6) displays all six archived
landmarks as grouped bars. For each landmark j, beta_j is divided by the 4.8.8
coordinate for that same landmark. Establishment uses the N18 reference and the
other five use the N15 response layer. No new coordinates, intervals or selections
are estimated. The short display labels beta_(up-arrow R) and beta_NSB identify
the archived interval-separated rebound coordinate and the coordinate of maximum
cumulative net strategy balance, respectively; NSB abbreviates net strategy
balance, (n_N_to_Q - n_Q_to_N)/N. The 4.8.8 ratio is exactly one by construction.
Horizontal offsets separate network categories for readability. The bars and
their endpoint markers do not represent temporal trajectories or an interpolated
parameter curve.

## 6. Single introduction and additional comparisons

Fig. S4 uses one periodic 4.8.8 graph with N = 33120, delta = 10, q_0 = 0.5,
global information and exactly one initial infection, rather than i_0 = 0.01.
There are 32 epidemics at beta_N = 2 and 4, and 128 at each of 6, 8, 12, 20, 32
and 50. All 832 reached absorption. Complete terminal outcomes were recovered;
832 dense temporal histories are not claimed.

| beta_N | Count with reach at least 5% |
|---|---|
| 2 | 0/32 |
| 4 | 0/32 |
| 6 | 9/128 |
| 8 | 68/128 |
| 12 | 93/128 |
| 20 | 110/128 |
| 32 | 86/128 |
| 50 | 30/128 |

At the 2% threshold, beta_N = 20 and 50 both have 110/128 successes; at 10%,
their counts are 99/128 and 8/128. The high-transmission decline depends on the
reach threshold. Fig. S4 A uses the 95% Wilson interval with
z = 1.959963984540054. For p = x/n, its center is
(p + z^2/(2n))/(1 + z^2/n), and half-width is
z sqrt(p(1-p)/n + z^2/(4n^2))/(1 + z^2/n). This is conditional on the fixed graph.
Panel B reclassifies the same epidemics at 2%, 5% and 10% reach.

Fig. S5 retains 96 adaptive/frozen matched pairs per beta_N on 4.8.8. Its spatial
bootstrap summary uses 5000 replicates. The RRG3 pilot uses four observed graphs
of N = 32768 and 16 pairs per graph at beta_N = 20 and 50, giving 256 epidemics
in total. The retained analysis in
`recovery/inputs/p8_original/code/analyze_p8.py` resamples four graphs, then 16
pairs per selected graph, for 10000 replicates using NumPy default_rng(20261012).
It preserves graph and replicate indices across the two beta values. Its stored
summary uses Adaptive-minus-Frozen; the figure uses Frozen-minus-Adaptive and
reverses endpoints accordingly. Both RRG3 contrast intervals include zero, with
mixed graph-level signs. Four graph clusters support an exploratory comparison.
The original construction seeds for those four graphs remain unavailable.

The larger-population comparisons in SI section 5, additional simulations in
section 6 and single-introduction experiments in section 7 remain separate.
The additional finite-population comparison used 12 realizations per setting,
rather than replacing the dense Fig. 4 estimates. Recorded histories were also
reclassified to examine detector sensitivity. These comparisons used the same
simulator and i_0 = 0.01. SI Table S3 preserves unresolved Square and heatmap
differences, incomplete RRG4 information checks, limited establishment and
high-order recurrence samples, and the unresolved additional Honeycomb/RRG3
contrast. New wording is not a resolution of those numerical differences.

## 7. Provenance and reconstruction scope

The figure index identifies measured input tables, archived source records and
the shared rendering modules. The source manifest covers the current versioned source tree; `DATA_BUNDLE.json`
and `ARTICLE_UPDATE.json` identify the base and current inputs. Historical
scientific observations and the simulation kernel are unchanged. Fig. 4 contains peak/duration and
Fig. 8 contains three landmarks. Supplementary Fig. S6 contains the relative
transmission bars; the active supplement has six figures. The axis windows in
section 9 remain unchanged. No epidemic or new bootstrap has been generated.

From the repository root, run `python main.py prepare-data --archive BASE.zip
--update-archive UPDATE.zip`, then `python main.py reproduce-published`.
The shared `fpm.visuals.article` module selects the current views through the
existing response/topology recipes. It filters/repositions axes and markers
without changing measurements or uncertainty endpoints. All18 active assets,
including the current S6 relative bars, are generated by this entry point.
Receipts and plotted-coordinate tables are supplied in `qa/` and
`article_update/qa/`. The full historical layouts remain preserved, outside the
active figure index. Current measured inputs are in `article_update/inputs/`.

The original network artwork in Fig. 1 is retained. ChatGPT and Codex assisted
code-based label standardization, assembly and other declared tasks; historical
model versions were not recorded. The scientific caption and manuscript retain
the applicable disclosure. Archived renderer names and execution labels provide
provenance only; they do not define separate physical dynamics.

Where a retained summary lacks a complete historical recipe or seed, that absence
must remain explicit. A replay of stored endpoints is not evidence of a newly
reconstructed historical calculation. Existing checks support their documented
scope and do not constitute new independent scientific acceptance.

## 8. Degree, contact organization and regime-dependent ordering

The new numerical comparisons in the narrative use the existing public data
archive with SHA-256
`96c4ee1193f35c3e95be1e92146f8f4faf120a04bb03b557d30ab050e6d2c8f0`.
They add interpretation and derived comparisons, not epidemic realizations or
new uncertainty intervals. All archive members below start with
`FPM_SHARED_20261102/`.

The general analysis in
[network_ordering_audit.py](../article_update/analysis/network_ordering_audit.py) reconstructs
reach as `ever_infected_count/population` from
`evidence/fig4/source/M27_G2_RUN_LEVEL_DERIVED.tsv`, averages within each graph
and then equally across graphs, and compares all 352 network reach means with
`reproduction/main/data/fig4_values.csv`. It covers 45,056 uncensored outcomes.
At delta = 10 and N = 32768, RRG4/RRG6 reach is 0.631053/0.698972 at beta_N = 4
and 0.734580/0.604254 at beta_N = 8. The manuscript rounds these means to three
decimal places. Archived marginal intervals are available; this analysis makes
no new simultaneous ordering test or estimate of a crossing location.

The 4,096 outcomes in `evidence/fig9/inputs/E3_RUNS_DERIVED.tsv` reconstruct the
arm means in `M28_E3_ARM_CELL_MEANS.csv` and the signed differences in
`reproduction/main/data/fig9_contrasts.csv`. The beta_N = 8 frozen RRG4/RRG6 means
are 0.822516/0.900839 and the adaptive means are 0.730568/0.567558. These controls
use 32 outcomes per arm on eight graphs, separately from the main 128-outcome
series. All eight RRG graphs per degree are reused; no run fingerprints overlap
with the main response series. Signed Frozen-minus-Adaptive ranks RRG6 > RRG4 >
RRG3 at beta_N = 8 and RRG3 > RRG4 > RRG6 at beta_N = 40. Its order is distinct
from the order of the individual-arm reaches or absolute contrast magnitudes.

The larger-population comparison comes from
`evidence/fig10/inputs/EXT3/05_G4_FINITE_SIZE/`:
`M32EXT3_G4_TRAJECTORY_CLASSIFICATION.csv` supplies reach, while `chunk0.tsv`
through `chunk3.tsv` supply the simulation design. At N = 131072, beta_N = 8 and
delta = 10, RRG4/RRG6 graph-first mean reach is 0.778421/0.777202 from 32
realizations and four graphs per degree. These are means, not the historical
summary medians. The analysis keeps this sample separate and tests no equivalence
or infinite-population persistence.

The landmark comparison reads `reproduction/main/data/fig8_landmarks.csv`.
No different-degree ordering violations occur within either family in its six
landmarks. The degree-three time-averaged-protection maxima at beta_N = 18
(4.8.8) and 19 (Honeycomb) have archived location intervals [18,22] and [17,19].
No new location interval is calculated. Establishment and later landmarks retain
their distinct populations and definitions.

From the repository root, reproduce these comparisons with:

```text
python article_update/analysis/network_ordering_audit.py --archive PATH_TO_PUBLIC_DATA_ZIP --output analysis/replay
```

The [retained audit receipt](../article_update/analysis/evidence/AUDIT_RECEIPT.json) identifies
the exact input members and hashes; its status refers to the preceding read-only
audit. The new replay has its own receipt. CSVs in `article_update/analysis/evidence/` retain the
values used for this revision without duplicating the large data archive. No
dynamic rule, event generator or per-figure simulation has been introduced.

## 9. Author-requested presentation windows, V5

The author requested Fig. 4 B through beta_N = 15, Fig. 8 through 20
and Fig. 11 through 15. Fig. 4 A remains through20.
The compositor removes drawn coordinates beyond the chosen maxima before adding
a small frame margin to keep boundary markers legible. It neither interpolates
new coordinates nor changes uncertainty endpoints. Fig. 8 retains all 24 selected
markers, its two linear segments and an internal legend. Fig. 11 retains five
pairs, y limits [-0.55, 0.30] and no error bars; all 220 source differences and
interval columns are retained in `article_update/inputs/fig11_contrasts.csv`, copied byte for
byte from the verified source. `article_update/qa/CURRENT_peak_duration.json`,
`CURRENT_selected_landmarks.json` and `CURRENT_degree_window.json` record the
actual visible coordinates and segments; their CSVs are in `article_update/qa/coordinates/`.
