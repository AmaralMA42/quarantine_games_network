# Core quantities and storage

The current manuscript and retained source code define the observables. These notes explain common columns; they do not assign invented physical units to unrecognized fields.

| Quantity or column | Meaning |
|---|---|
| `population`, N | Number of vertices or individuals in the simulation. |
| `beta_n`, beta_N | Normal-strategy transmission parameter. |
| `beta_q`, beta_Q | Protected-strategy transmission parameter. |
| `beta_a` | Cross-strategy transmission parameter resolved by the recorded mixing contract. |
| `beta_mix_a` | Mixing coefficient in that contract. |
| `delta` | Perceived-risk-strength parameter. |
| `tau` | Strategy-update rate; historical `chi` is normalized at ingestion by the retained contract. |
| `fermi_k` | Fermi comparison softness; historical `epsilon` is an ingestion alias. |
| `i0`, `q0` | Initial infected and protected fractions. Individual-introduction records use I0 = 1. |
| `ever_infected_count` | Number ever infected, including initially infected individuals. |
| `reach`, R_infinity | `ever_infected_count / population`. |
| `peak_infected` | Maximum infected count in the retained record. |
| `imax_fraction`, I_max | Peak infected count divided by population. |
| `duration`, T_epi | Epidemic duration to the specified stopping condition, in model-time units. |
| `q_final_fraction`, q_f | Final active protected fraction, excluding recovered individuals. |
| `q_time_average` | Mean active protection over epidemic duration, as defined by the retained protection-time integral. |
| `graph_sha256` | Identity of the observed graph artifact; not a numerical observable. |
| `run_fingerprint` | Identity of the recorded scientific run specification. |
| `*_seed` | Recorded random-stream metadata. Read as strings or exact integers; do not coerce 64-bit identifiers through floating point. |

The payoff-equality prevalence I_star is distinct from the peak infected fraction I_max. Peak-count probabilities and conditional uncertainty follow the detector and sampling definitions retained with each dataset; the detector used for measured heatmaps is not silently replaced by the Figure 10 detector. Lines in parameter-response figures connect measured coordinates. Empty matrix coordinates are not measured zeros.

The column dictionary describes sampled storage forms. It does not certify full precision of every historical seed field, infer units from numeric values, or turn all run-level fields into independent samples. Graph-first dependence and shared RRG banks remain relevant when combining network contrasts.
