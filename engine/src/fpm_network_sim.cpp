// FPM_NETWORK_SIM_V5_FERMI_K_NOTATION_CANDIDATE
// Generic Quarantine Games simulator. One executable consumes graph artifacts and run specs;
// topology and scientific stage are data, never separate simulator programs.
// Scientific base: M20B corrected independent fixed-count Q/I initialization.
// Engineering changes are regression-gated against the M20B as-run kernel.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <charconv>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <map>


#include "fpm/stable_rng.hpp"
#include "fpm/pools.hpp"
#include "fpm/periodic_certificate.hpp"
#include "fpm/graph.hpp"
#include "fpm/run_spec.hpp"
#include "fpm/beta_contract.hpp"
#include "fpm/outcome.hpp"
#include "fpm/sha256.hpp"

using fpm::StableRng;
using fpm::Pool;
using fpm::Fenwick;
using fpm::EdgePool;
using fpm::EdgePoolBackend;
using fpm::Graph;
using fpm::GraphRuntime;
using fpm::build_graph_runtime;
using fpm::Config;
using fpm::Outcome;

// Epidemiological/behavioural product states.
// N/Q encode the social strategy while S/I/R encode disease status.  Recovery
// is absorbing epidemiologically, so one RR state is sufficient for the current
// Quarantine Games branch.  These integer values are serialized indirectly via
// observables and must not be reordered without a schema/version change.
enum State : int { SN = 0, SQ = 1, IN = 2, IQ = 3, RR = 4 };



static bool is_lower_hex_sha256(const std::string &text) {
    if (text.size() != 64) return false;
    return std::all_of(text.begin(), text.end(), [](unsigned char ch) {
        return (ch >= '0' && ch <= '9') || (ch >= 'a' && ch <= 'f');
    });
}

static std::string json_escape_ascii(const std::string &text) {
    std::ostringstream out; out << '"';
    const char *hex = "0123456789abcdef";
    for (unsigned char ch : text) {
        if (ch > 0x7f) throw std::runtime_error("semantic identity strings must be ASCII for FPM_TSV_ID_V3");
        switch (ch) {
            case '"': out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\b': out << "\\b"; break;
            case '\f': out << "\\f"; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (ch < 0x20) out << "\\u00" << hex[(ch >> 4) & 0xf] << hex[ch & 0xf];
                else out << char(ch);
        }
    }
    out << '"'; return out.str();
}

static std::string python_float_json(const std::string &raw) {
    const double value = std::stod(raw);
    if (!std::isfinite(value)) throw std::runtime_error("non-finite semantic identity float");
    std::array<char,128> buffer{};
    const auto result = std::to_chars(buffer.data(), buffer.data()+buffer.size(), value, std::chars_format::general);
    if (result.ec != std::errc{}) throw std::runtime_error("could not canonicalize semantic identity float");
    std::string text(buffer.data(), result.ptr);
    if (text.find_first_of(".eE") == std::string::npos) text += ".0";
    return text;
}

static std::string canonical_semantic_identity_json(
    const std::vector<std::string> &fields,
    const std::unordered_map<std::string,int> &column
) {
    static const std::unordered_set<std::string> excluded{
        "run_fingerprint", "run_index", "run_id", "graph_path",
        "graph_certificate_path", "fingerprint_version"
    };
    static const std::unordered_set<std::string> integer_fields{
        "population","n0","n_occ","n_gc","graph_seed","q_init_seed","i_init_seed",
        "dynamic_seed","signal_seed","information_radius","event_limit"
    };
    static const std::unordered_set<std::string> unsigned_integer_fields{
        "graph_seed","q_init_seed","i_init_seed","dynamic_seed","signal_seed"
    };
    static const std::unordered_set<std::string> float_fields{
        "beta_n","beta_a","beta_q","beta_mix_a","q0","i0","tau","chi","i_c","fermi_k","epsilon","omega","delta",
        "information_alpha","behavior_bin_width","tmax"
    };
    static const std::vector<std::string> required{
        "schema_version","analysis_role","run_fingerprint","run_index","run_id","stage",
        "topology","graph_path","graph_sha256","graph_certificate_path",
        "graph_certificate_sha256","boundary_condition","population","n0","n_occ","n_gc",
        "graph_seed","q_init_seed","i_init_seed","dynamic_seed","contact_mode","beta_n","beta_a","beta_q","beta_mix_a","beta_mixing_rule",
        "profile","arm","protocol","q0","i0","tau","i_c","fermi_k",
        "omega","delta","seed_policy","behavior_mode","q_placement","information_mode",
        "information_radius","information_alpha","signal_seed","instrumentation_mode",
        "trace_mode","behavior_bin_width","event_limit","tmax"
    };
    std::map<std::string,std::string> entries;
    for (const auto &canonical_key : required) {
        std::string key = canonical_key;
        if (canonical_key == "tau" && column.find("tau") == column.end()
            && column.find("chi") != column.end()) key = "chi";
        if (canonical_key == "fermi_k" && column.find("fermi_k") == column.end()
            && column.find("epsilon") != column.end()) key = "epsilon";
        if (excluded.count(key)) continue;
        const auto found=column.find(key);
        if (found==column.end()) throw std::runtime_error("semantic identity missing field: "+key);
        const auto &raw=fields.at(found->second);
        if (integer_fields.count(key)) {
            entries[key] = unsigned_integer_fields.count(key)
                ? std::to_string(std::stoull(raw)) : std::to_string(std::stoll(raw));
        } else if (float_fields.count(key)) {
            entries[key] = python_float_json(raw);
        } else {
            entries[key] = json_escape_ascii(raw);
        }
    }
    std::ostringstream out; out << '{'; bool first=true;
    for (const auto &[key,value] : entries) {
        if (!first) out << ',';
        first=false;
        out << json_escape_ascii(key) << ':' << value;
    }
    out << '}'; return out.str();
}

static std::vector<std::string> split(const std::string &line, char delimiter) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string item;
    while (std::getline(stream, item, delimiter)) fields.push_back(item);
    return fields;
}

// Bit-exact SplitMix64 stream plus rejection sampling.  Unlike std::shuffle
// and std::*_distribution, this mapping is part of the Foundation contract.
// Fermi imitation probability.  Saturating before exp() avoids overflow while
// preserving the historical numerical contract in the tails.
static double fermi(double difference, double fermi_k) {
    const double x = difference / fermi_k;
    if (x >= 40.0) return 1.0;
    if (x <= -40.0) return 0.0;
    return 1.0 / (1.0 + std::exp(-x));
}

// Build the radius-r information neighborhood used by local-information modes.
// This is graph-invariant and should be shared/reused whenever those modes are
// promoted to the stable surface.  The current corrected M20B surface is global.
static std::vector<std::vector<int>> build_signal_neighborhood(
    const Graph &graph, int radius
) {
    if (radius == 1) return graph.neighbors;
    std::vector<std::vector<int>> neighborhoods(graph.n);
    std::vector<int> seen(graph.n, -1);
    std::vector<int> distance(graph.n, -1);
    std::vector<int> queue;
    queue.reserve(graph.n);
    for (int root = 0; root < graph.n; ++root) {
        queue.clear();
        queue.push_back(root);
        seen[root] = root;
        distance[root] = 0;
        size_t head = 0;
        while (head < queue.size()) {
            const int u = queue[head++];
            if (distance[u] == radius) continue;
            for (int v : graph.neighbors[u]) {
                if (seen[v] == root) continue;
                seen[v] = root;
                distance[v] = distance[u] + 1;
                queue.push_back(v);
                neighborhoods[root].push_back(v);
            }
        }
    }
    return neighborhoods;
}

// Known regular production topologies.  The graph file remains authoritative;
// this helper is a validation contract, not a topology-specific dynamics switch.
static int expected_degree(const std::string &topology) {
    if (topology == "honeycomb" || topology == "rr_k3" || topology == "truncated_square_488" || topology == "4.8.8") return 3;
    if (topology == "square" || topology == "kagome" || topology == "rr_k4") return 4;
    if (topology == "triangular" || topology == "rr_k6") return 6;
    throw std::runtime_error("unknown topology: " + topology);
}

// Parse and validate a simple undirected graph from the canonical text format.
// The simulator never generates scientific graph realizations implicitly: graph
// materiality and hashes belong to the run specification/provenance layer.
static Graph load_graph(const std::string &path, const std::string &topology) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("could not open graph: " + path);
    Graph graph;
    int m = 0;
    input >> graph.n >> m;
    if (graph.n < 2 || m < 1) throw std::runtime_error("invalid graph header");
    graph.neighbors.resize(graph.n);
    graph.degree.assign(graph.n, 0);
    graph.edges.reserve(m);
    std::unordered_set<uint64_t> unique_edges;
    unique_edges.reserve(size_t(m) * 2);
    for (int e = 0; e < m; ++e) {
        int u, v;
        input >> u >> v;
        if (!input || u < 0 || v < 0 || u >= graph.n || v >= graph.n || u == v) {
            throw std::runtime_error("invalid graph edge");
        }
        const uint32_t a = uint32_t(std::min(u, v));
        const uint32_t b = uint32_t(std::max(u, v));
        const uint64_t key = (uint64_t(a) << 32) | uint64_t(b);
        if (!unique_edges.insert(key).second) {
            throw std::runtime_error("duplicate graph edge");
        }
        graph.edges.push_back({u, v});
        graph.neighbors[u].push_back(v);
        graph.neighbors[v].push_back(u);
        graph.degree[u]++;
        graph.degree[v]++;
    }
    std::string trailing;
    if (input >> trailing) throw std::runtime_error("trailing graph records");
    const int degree_expected = expected_degree(topology);
    for (int degree : graph.degree) {
        if (degree != degree_expected) throw std::runtime_error("graph degree mismatch");
    }
    std::vector<char> seen(graph.n, 0);
    std::vector<int> queue{0};
    seen[0] = 1;
    for (size_t head = 0; head < queue.size(); ++head) {
        for (int neighbor : graph.neighbors[queue[head]]) {
            if (!seen[neighbor]) {
                seen[neighbor] = 1;
                queue.push_back(neighbor);
            }
        }
    }
    if (int(queue.size()) != graph.n) throw std::runtime_error("disconnected graph");
    std::vector<int> order(graph.n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return graph.degree[a] > graph.degree[b];
    });
    graph.hub.assign(graph.n, 0);
    const int hub_count = std::max(1, int(std::ceil(0.01 * graph.n)));
    for (int i = 0; i < hub_count; ++i) graph.hub[order[i]] = 1;
    return graph;
}

// Select one CTMC event family from cumulative hazards.  The exact ordering is
// part of the historical seed-to-trajectory contract for the fenwick_exact path.
static int choose_event(const std::array<double, 9> &rates, double draw) {
    double cumulative = 0.0;
    for (int i = 0; i < int(rates.size()); ++i) {
        cumulative += rates[i];
        if (draw < cumulative) return i;
    }
    return int(rates.size()) - 1;
}

// Convert perceived prevalence into N->Q and Q->N imitation probabilities.
// This is behaviour, not biological transmission; keep the two mechanisms
// separate in both code and interpretation.
static std::pair<double, double> adoption_from_prevalence(
    const Config &cfg, double prevalence
) {
    prevalence = std::clamp(prevalence, 0.0, 1.0);
    const double pi_q = -cfg.omega;
    const double beta_n = cfg.beta_n;
    const double pi_n = -cfg.delta * beta_n * prevalence;
    return {
        fermi(pi_q - pi_n, cfg.fermi_k),
        fermi(pi_n - pi_q, cfg.fermi_k)
    };
}

static std::pair<double, double> adoption(const Config &cfg, int infected) {
    return adoption_from_prevalence(
        cfg, double(infected) / double(cfg.population)
    );
}

// Well-mixed reference process sharing the corrected initialization/RNG contract.
// It is a formalism branch inside the same executable, not a substitute for the
// network model when geometry is part of the scientific question.
static Outcome run_well_mixed(const Config &cfg) {
    const int n = cfg.population;
    StableRng i_rng(cfg.i_init_seed);
    StableRng rng(cfg.dynamic_seed);
    const int initial_i = std::clamp(
        cfg.i0 > 0.0 ? int(std::llround(cfg.i0 * n)) : 1, 1, n
    );
    const int q_total = std::clamp(int(std::llround(cfg.q0 * n)), 0, n);
    int initial_iq = 0;
    int q_remaining = q_total;
    int population_remaining = n;
    for (int draw = 0; draw < initial_i; ++draw) {
        if (i_rng.bounded(uint64_t(population_remaining)) < uint64_t(q_remaining)) {
            initial_iq++;
            q_remaining--;
        }
        population_remaining--;
    }
    const int initial_in = initial_i - initial_iq;
    const int initial_sq = std::clamp(q_total - initial_iq, 0, n - initial_i);
    const int initial_sn = n - initial_i - initial_sq;
    std::array<int64_t, 5> count{initial_sn, initial_sq, initial_in, initial_iq, 0};
    if (count[SN] < 0 || count[SQ] < 0) {
        throw std::runtime_error("bad PAPER_EQUIVALENT initialization");
    }
    Outcome out;
    out.initial_i = initial_i;
    out.initial_q = q_total;
    out.initial_iq = initial_iq;
    out.ever_infected_count = initial_i;
    out.seed_vertex = int(cfg.dynamic_seed % uint64_t(n));
    out.peak_infected = int(count[IN] + count[IQ]);
    out.max_q = int(count[SQ] + count[IQ]);
    out.q_mesoscopic = out.max_q >= int(std::ceil(std::sqrt(double(n))));
    out.q_macroscopic = double(out.max_q) / n >= 0.05;
    out.activated_ic = double(out.peak_infected) / n >= cfg.i_c;
    const double beta_n = cfg.beta_n;
    const double beta_q = cfg.beta_q;
    const double beta_a = cfg.beta_a;
    const double denom = double(n - 1);
    double time = 0.0;
    while (count[IN] + count[IQ] > 0 && out.events < cfg.event_limit && time < cfg.tmax) {
        const int infected = int(count[IN] + count[IQ]);
        const auto p = cfg.behavior_mode == "adaptive"
            ? adoption(cfg, infected)
            : std::pair<double, double>{0.0, 0.0};
        const double active_n = double(count[SN] + count[IN]);
        const double active_q = double(count[SQ] + count[IQ]);
        std::array<double, 9> rates{
            beta_n * count[SN] * count[IN] / denom,
            beta_a * count[SN] * count[IQ] / denom,
            beta_a * count[SQ] * count[IN] / denom,
            beta_q * count[SQ] * count[IQ] / denom,
            cfg.tau * count[SN] * active_q * p.first / denom,
            cfg.tau * count[SQ] * active_n * p.second / denom,
            cfg.tau * count[IN] * active_q * p.first / denom,
            cfg.tau * count[IQ] * active_n * p.second / denom,
            double(infected)
        };
        const double total = std::accumulate(rates.begin(), rates.end(), 0.0);
        if (!(total > 0.0) || !std::isfinite(total)) break;
        const double previous_time = time;
        const double next_time = time
            + -std::log(std::max(1e-300, rng.uniform01())) / total;
        time = std::min(next_time, cfg.tmax);
        out.q_person_time += double(count[SQ] + count[IQ]) / n
            * (std::min(time, cfg.tmax) - previous_time);
        if (next_time > cfg.tmax) break;
        const int previous_q = int(count[SQ] + count[IQ]);
        const int event = choose_event(rates, rng.uniform01() * total);
        switch (event) {
            case 0: case 1: count[SN]--; count[IN]++; out.ever_infected_count++; break;
            case 2: case 3: count[SQ]--; count[IQ]++; out.ever_infected_count++; break;
            case 4: count[SN]--; count[SQ]++; break;
            case 5: count[SQ]--; count[SN]++; break;
            case 6: count[IN]--; count[IQ]++; break;
            case 7: count[IQ]--; count[IN]++; break;
            case 8:
                if (count[IQ] == 0 || (count[IN] > 0 && rng.uniform01() * infected < count[IN])) {
                    count[IN]--; count[RR]++;
                } else {
                    count[IQ]--; count[RR]++;
                }
                break;
        }
        out.events++;
        out.accepted_event_count++;
        const int current_i = int(count[IN] + count[IQ]);
        if (current_i > out.peak_infected) {
            out.peak_infected = current_i;
            out.peak_time = time;
        }
        out.max_q = std::max(out.max_q, int(count[SQ] + count[IQ]));
        if (previous_q > 0 && count[SQ] + count[IQ] == 0) out.q_extinction_time = time;
        out.q_mesoscopic = out.q_mesoscopic
            || out.max_q >= int(std::ceil(std::sqrt(double(n))));
        out.q_macroscopic = out.q_macroscopic || double(out.max_q) / n >= 0.05;
        out.activated_ic = out.activated_ic || double(current_i) / n >= cfg.i_c;
    }
    out.censored = (count[IN] + count[IQ] > 0);
    out.outbreak_size = out.ever_infected_count;
    out.duration = time;
    out.stop_reason = out.censored
        ? (time >= cfg.tmax ? "tmax" : "event_limit")
        : "absorbed_I0";
    out.final_q = int(count[SQ] + count[IQ]);
    if (out.max_q > 0 && out.final_q == 0 && !std::isfinite(out.q_extinction_time)) {
        out.q_extinction_time = time;
    }
    return out;
}

// One stochastic network trajectory over one immutable GraphRuntime.
//
// Design invariants:
//   * graph-only indexing lives in GraphRuntime and is reused across trajectories;
//   * q_rng, i_rng and rng are domain-separated streams;
//   * fenwick_exact is the reference backend and preserves historical mapping;
//   * uniform_regular is opt-in and may preserve the law without preserving the
//     same edge choice for a given seed, so it is not claim-grade by default;
//   * observers may read state but must not consume the dynamics RNG stream.
class GraphProcess {
    const Config &cfg;
    const Graph &graph;
    const GraphRuntime &runtime;
    EdgePoolBackend edge_backend;
    std::vector<int> state;
    std::array<int64_t, 5> count{};
    // Infection and social-imitation hazards are independent competing clocks.
    // A directed edge may therefore belong to one pool in each family at once.
    std::vector<int> infection_edge_category;
    std::vector<int> behavior_edge_category;
    std::array<EdgePool, 8> pools;
    std::array<int64_t, 8> raw_pool_counts{};
    std::array<Fenwick, 4> behavior_node_pools;
    std::vector<int> behavior_node_category;
    std::vector<double> behavior_node_weight;
    std::vector<int> adjacent_infected_count;
    std::vector<int> signal_infected_count;
    const std::vector<std::vector<int>> *signal_neighborhood = nullptr;
    std::vector<int> signal_donor;
    std::vector<int> signal_recipient;
    std::vector<int> distance_from_seed;
    std::vector<int> shell_total_count;
    std::vector<int> shell_q_count;
    Pool infected;
    std::vector<char> ever_infected;
    int ever_infected_count = 0;
    int64_t infected_degree_sum = 0;
    int hubs_infected = 0;
    StableRng q_rng;
    StableRng i_rng;
    StableRng rng;
    std::unordered_map<int, int> qn_bins;
    std::unordered_map<int, int> nq_bins;
    std::unordered_map<int, int> i_bin_max;
    std::unordered_map<int, int> q_bin_last;
    struct TraceSnapshot {
        std::array<int64_t, 5> state_count{};
        int64_t si_n_edges = 0;
        int64_t si_q_edges = 0;
        int q_interface_nodes = 0;
        int q_far_nodes = 0;
        double infected_radius_mean = 0.0;
        double infected_radius_rms = 0.0;
        int max_ever_distance = 0;
        std::vector<int> shell_q_count;
    };
    std::unordered_map<int, TraceSnapshot> detail_bins;
    int q_interface_nodes = 0;
    int cached_global_infected = -1;
    std::pair<double, double> cached_global_adoption{0.0, 0.0};
    std::vector<int> interface_affected_scratch;
    std::vector<int> local_affected_scratch;
    int max_ever_distance = 0;
    int ever_distance_count = 0;
    double ever_distance_square_sum = 0.0;
    double infected_distance_sum = 0.0;
    double infected_distance_square_sum = 0.0;
    struct TraceV2Point { double time; int infected; int q; };
    std::vector<TraceV2Point> trace_v2_points;
    double next_trace_time = 0.0;

    static bool is_infected(int s) { return s == IN || s == IQ; }
    static bool is_q(int s) { return s == SQ || s == IQ; }
    static bool is_n(int s) { return s == SN || s == IN; }

    const std::pair<double, double> &global_adoption_cached(int infected_count) {
        if (infected_count != cached_global_infected) {
            cached_global_infected = infected_count;
            cached_global_adoption = adoption(cfg, infected_count);
        }
        return cached_global_adoption;
    }

    int infection_category(int u, int v) const {
        const int a = state[u], b = state[v];
        if (a == SN && b == IN) return 0;
        if (a == SN && b == IQ) return 1;
        if (a == SQ && b == IN) return 2;
        if (a == SQ && b == IQ) return 3;
        return -1;
    }

    int behavior_category(int u, int v) const {
        const int a = state[u], b = state[v];
        if (a == SN && is_q(b)) return 4;
        if (a == SQ && is_n(b)) return 5;
        if (a == IN && is_q(b)) return 6;
        if (a == IQ && is_n(b)) return 7;
        return -1;
    }

    double weight_for(int cat, int u) const {
        if (cat < 0) return 0.0;
        if (cat < 4 && cfg.contact_mode == "per_edge") return 1.0;
        return 1.0 / double(graph.degree[u]);
    }

    void refresh(int edge) {
        const int u = runtime.directed_u[edge];

        const int old_infection = infection_edge_category[edge];
        if (old_infection >= 0) {
            pools[old_infection].add(edge, -weight_for(old_infection, u));
            raw_pool_counts[old_infection]--;
        }
        const int old_behavior = behavior_edge_category[edge];
        if (old_behavior >= 0) {
            pools[old_behavior].add(edge, -weight_for(old_behavior, u));
            raw_pool_counts[old_behavior]--;
        }

        const int next_infection = infection_category(u, runtime.directed_v[edge]);
        int next_behavior = behavior_category(u, runtime.directed_v[edge]);
        // Exact historical reproduction keeps the legacy single-category edge
        // semantics. Corrected claim-bearing science uses independent competing
        // infection and imitation clocks on the same directed contact.
        if (cfg.beta_mixing_rule == "LEGACY_PUBLISHED" && next_infection >= 0) {
            next_behavior = -1;
        }
        infection_edge_category[edge] = next_infection;
        behavior_edge_category[edge] = next_behavior;

        if (next_infection >= 0) {
            pools[next_infection].add(edge, weight_for(next_infection, u));
            raw_pool_counts[next_infection]++;
        }
        if (next_behavior >= 0) {
            pools[next_behavior].add(edge, weight_for(next_behavior, u));
            raw_pool_counts[next_behavior]++;
        }
    }

    bool has_local_signal() const {
        return cfg.information_mode == "local"
            || cfg.information_mode == "local_shuffled"
            || cfg.information_mode == "mixed";
    }

    bool uses_node_behavior_pools() const {
        return cfg.information_mode == "local"
            || cfg.information_mode == "local_shuffled";
    }

    bool uses_mixed_thinning() const {
        return cfg.information_mode == "mixed";
    }

    static int node_behavior_category(int s) {
        if (s == SN) return 0;
        if (s == SQ) return 1;
        if (s == IN) return 2;
        if (s == IQ) return 3;
        return -1;
    }

    double raw_local_prevalence(int u) const {
        return signal_neighborhood != nullptr
            && !(*signal_neighborhood)[u].empty()
            ? double(signal_infected_count[u])
                / double((*signal_neighborhood)[u].size())
            : 0.0;
    }

    double information_signal(int u) const {
        const int donor = signal_donor.empty() ? u : signal_donor[u];
        const double local = raw_local_prevalence(donor);
        if (cfg.information_mode != "mixed") return local;
        const double global = double(infected.size()) / double(graph.n);
        return cfg.information_alpha * global
            + (1.0 - cfg.information_alpha) * local;
    }

    double focal_local_prevalence(int u) const {
        if (graph.degree[u] <= 0) return 0.0;
        int local_infected = 0;
        for (int v : graph.neighbors[u]) local_infected += is_infected(state[v]);
        return double(local_infected) / double(graph.degree[u]);
    }

    void refresh_behavior_node(int u) {
        if (!uses_node_behavior_pools()) return;
        const int old_category = behavior_node_category[u];
        if (old_category >= 0) {
            behavior_node_pools[old_category].add(u, -behavior_node_weight[u]);
        }
        int next_category = -1;
        double next_weight = 0.0;
        if (cfg.behavior_mode == "adaptive" && state[u] != RR) {
            next_category = node_behavior_category(state[u]);
            if (graph.degree[u] > 0 && is_n(state[u])) {
                int opposite_q = 0;
                for (int v : graph.neighbors[u]) opposite_q += is_q(state[v]);
                const auto focal_p =
                    adoption_from_prevalence(cfg, information_signal(u));
                next_weight = cfg.tau * focal_p.first
                    * double(opposite_q) / double(graph.degree[u]);
            } else if (graph.degree[u] > 0 && is_q(state[u])) {
                double model_weight = 0.0;
                for (int v : graph.neighbors[u]) {
                    if (!is_n(state[v])) continue;
                    const auto model_p =
                        adoption_from_prevalence(cfg, information_signal(v));
                    model_weight += model_p.second;
                }
                next_weight = cfg.tau * model_weight / double(graph.degree[u]);
            }
        }
        behavior_node_category[u] = next_category;
        behavior_node_weight[u] = next_weight;
        if (next_category >= 0) {
            behavior_node_pools[next_category].add(u, next_weight);
        }
    }

    void initialize_information_and_interface() {
        adjacent_infected_count.assign(graph.n, 0);
        for (int u = 0; u < graph.n; ++u) {
            for (int v : graph.neighbors[u]) {
                adjacent_infected_count[u] += is_infected(state[v]);
            }
        }
        q_interface_nodes = 0;
        for (int u = 0; u < graph.n; ++u) {
            q_interface_nodes += is_q(state[u]) && adjacent_infected_count[u] > 0;
        }
        if (!has_local_signal()) return;
        if (signal_neighborhood == nullptr) {
            throw std::runtime_error("missing local signal neighborhood");
        }
        signal_infected_count.assign(graph.n, 0);
        for (int u = 0; u < graph.n; ++u) {
            for (int v : (*signal_neighborhood)[u]) {
                signal_infected_count[u] += is_infected(state[v]);
            }
        }
        signal_donor.resize(graph.n);
        signal_recipient.resize(graph.n);
        std::iota(signal_donor.begin(), signal_donor.end(), 0);
        if (cfg.information_mode == "local_shuffled") {
            StableRng signal_rng(cfg.signal_seed);
            signal_rng.shuffle(signal_donor);
        }
        for (int recipient = 0; recipient < graph.n; ++recipient) {
            signal_recipient[signal_donor[recipient]] = recipient;
        }
        if (uses_node_behavior_pools()) {
            for (int u = 0; u < graph.n; ++u) refresh_behavior_node(u);
        }
    }

    void mark_ever(int u) {
        if (ever_infected[u]) return;
        ever_infected[u] = 1;
        ever_infected_count++;
        infected_degree_sum += graph.degree[u];
        hubs_infected += graph.hub[u];
        if (!distance_from_seed.empty() && distance_from_seed[u] >= 0) {
            const int distance = distance_from_seed[u];
            max_ever_distance = std::max(max_ever_distance, distance);
            ever_distance_count++;
            ever_distance_square_sum += double(distance) * distance;
        }
    }

    void initialize_trace_v2() {
        if (!cfg.trace_mode) return;
        trace_v2_points.push_back({0.0, int(infected.size()), int(count[SQ] + count[IQ])});
        next_trace_time = cfg.behavior_bin_width;
    }

    void fill_trace_v2(double horizon, int current_i, int current_q) {
        if (!cfg.trace_mode) return;
        const double tolerance = 1e-12 * std::max(1.0, horizon);
        while (next_trace_time <= horizon + tolerance) {
            trace_v2_points.push_back({next_trace_time, current_i, current_q});
            next_trace_time += cfg.behavior_bin_width;
        }
    }

    std::string serialize_trace_v2(double final_time, int final_i, int final_q) const {
        if (!cfg.trace_mode) return "";
        std::ostringstream stream;
        stream << std::setprecision(17);
        for (size_t index = 0; index < trace_v2_points.size(); ++index) {
            if (index) stream << ';';
            const auto &point = trace_v2_points[index];
            stream << point.time << ':' << point.infected << ':' << point.q;
        }
        stream << ";F:" << final_time << ':' << final_i << ':' << final_q;
        return stream.str();
    }

    void change_state(int u, int next) {
        const int previous = state[u];
        if (previous == next) return;
        const bool infection_changed = is_infected(previous) != is_infected(next);
        const bool q_changed = is_q(previous) != is_q(next);
        interface_affected_scratch.clear();
        interface_affected_scratch.reserve(graph.degree[u] + 1);
        interface_affected_scratch.push_back(u);
        for (int v : graph.neighbors[u]) interface_affected_scratch.push_back(v);
        for (int v : interface_affected_scratch) {
            q_interface_nodes -= is_q(state[v]) && adjacent_infected_count[v] > 0;
        }
        count[previous]--;
        count[next]++;
        if (is_infected(previous)) {
            infected.remove(u);
            if (!distance_from_seed.empty() && distance_from_seed[u] >= 0) {
                infected_distance_sum -= distance_from_seed[u];
                infected_distance_square_sum -=
                    double(distance_from_seed[u]) * distance_from_seed[u];
            }
        }
        state[u] = next;
        if (q_changed && !distance_from_seed.empty() && distance_from_seed[u] >= 0) {
            const int shell = distance_from_seed[u];
            if (shell < int(shell_q_count.size())) {
                shell_q_count[shell] += is_q(next) ? 1 : -1;
            }
        }
        if (is_infected(next)) {
            infected.add(u);
            if (!distance_from_seed.empty() && distance_from_seed[u] >= 0) {
                infected_distance_sum += distance_from_seed[u];
                infected_distance_square_sum +=
                    double(distance_from_seed[u]) * distance_from_seed[u];
            }
            mark_ever(u);
        }
        for (int pos = runtime.incident_offsets[u]; pos < runtime.incident_offsets[u + 1]; ++pos) {
            const int e = runtime.incident_edges[pos];
            refresh(e);
            refresh(e ^ 1);
        }
        if (infection_changed) {
            const int delta = is_infected(next) ? 1 : -1;
            for (int v : graph.neighbors[u]) adjacent_infected_count[v] += delta;
        }
        for (int v : interface_affected_scratch) {
            q_interface_nodes += is_q(state[v]) && adjacent_infected_count[v] > 0;
        }
        if (has_local_signal()) {
            local_affected_scratch.clear();
            local_affected_scratch.reserve(
                graph.degree[u] + (has_local_signal()
                    ? (*signal_neighborhood)[u].size() : 0) + 2
            );
            local_affected_scratch.push_back(u);
            for (int v : graph.neighbors[u]) local_affected_scratch.push_back(v);
            if (infection_changed) {
                const int delta = is_infected(next) ? 1 : -1;
                for (int focal : (*signal_neighborhood)[u]) {
                    signal_infected_count[focal] += delta;
                    const int recipient =
                        cfg.information_mode == "local_shuffled"
                            ? signal_recipient[focal] : focal;
                    local_affected_scratch.push_back(recipient);
                    // Under strict pairwise LOCAL semantics, a changed N-model payoff
                    // also changes every neighboring Q focal's Q->N imitation hazard.
                    for (int neighbor : graph.neighbors[recipient]) {
                        local_affected_scratch.push_back(neighbor);
                    }
                }
            }
            std::sort(local_affected_scratch.begin(), local_affected_scratch.end());
            local_affected_scratch.erase(
                std::unique(local_affected_scratch.begin(), local_affected_scratch.end()),
                local_affected_scratch.end()
            );
            if (uses_node_behavior_pools()) {
                for (int v : local_affected_scratch) refresh_behavior_node(v);
            }
        }
    }

    int select_seed() {
        std::vector<int> candidates;
        if (cfg.seed_policy == "uniform") {
            return int(i_rng.bounded(uint64_t(graph.n)));
        }
        std::vector<int> degrees = graph.degree;
        std::sort(degrees.begin(), degrees.end());
        if (cfg.seed_policy == "low") {
            const int cutoff = degrees[(graph.n - 1) / 2];
            for (int u = 0; u < graph.n; ++u) if (graph.degree[u] <= cutoff) candidates.push_back(u);
        } else if (cfg.seed_policy == "typical") {
            const int low_cutoff = degrees[int(std::floor(0.40 * (graph.n - 1)))];
            const int high_cutoff = degrees[int(std::ceil(0.60 * (graph.n - 1)))];
            for (int u = 0; u < graph.n; ++u) {
                if (graph.degree[u] >= low_cutoff && graph.degree[u] <= high_cutoff) {
                    candidates.push_back(u);
                }
            }
        } else if (cfg.seed_policy == "hub") {
            for (int u = 0; u < graph.n; ++u) if (graph.hub[u]) candidates.push_back(u);
        } else {
            throw std::runtime_error("unknown seed policy");
        }
        return candidates[size_t(i_rng.bounded(candidates.size()))];
    }

    int select_q_seed(int infection_seed) {
        if (cfg.q_placement == "none") return -1;
        if (cfg.q_placement == "adjacent") {
            const auto &neighbors = graph.neighbors[infection_seed];
            if (neighbors.empty()) throw std::runtime_error("infection seed has no neighbor");
            return neighbors[size_t(q_rng.bounded(neighbors.size()))];
        }
        if (cfg.q_placement == "uniform") {
            int candidate = infection_seed;
            while (candidate == infection_seed) {
                candidate = int(q_rng.bounded(uint64_t(graph.n)));
            }
            return candidate;
        }
        throw std::runtime_error("unknown Q placement");
    }

    void prepare_seed_distances(int seed) {
        distance_from_seed.assign(graph.n, -1);
        std::vector<int> queue;
        queue.reserve(graph.n);
        queue.push_back(seed);
        distance_from_seed[seed] = 0;
        size_t head = 0;
        while (head < queue.size()) {
            const int u = queue[head++];
            for (int v : graph.neighbors[u]) {
                if (distance_from_seed[v] >= 0) continue;
                distance_from_seed[v] = distance_from_seed[u] + 1;
                queue.push_back(v);
            }
        }
        const int max_shell = *std::max_element(
            distance_from_seed.begin(), distance_from_seed.end()
        );
        shell_total_count.assign(max_shell + 1, 0);
        shell_q_count.assign(max_shell + 1, 0);
        for (int distance : distance_from_seed) {
            if (distance >= 0) shell_total_count[distance]++;
        }
    }

    void capture_detail_bin(int bin) {
        if (!cfg.trace_mode) return;
        TraceSnapshot snapshot;
        snapshot.state_count = count;
        snapshot.si_n_edges = raw_pool_counts[0] + raw_pool_counts[1];
        snapshot.si_q_edges = raw_pool_counts[2] + raw_pool_counts[3];
        snapshot.q_interface_nodes = q_interface_nodes;
        snapshot.q_far_nodes =
            int(count[SQ] + count[IQ]) - q_interface_nodes;
        if (infected.size() > 0) {
            snapshot.infected_radius_mean =
                infected_distance_sum / double(infected.size());
            snapshot.infected_radius_rms = std::sqrt(
                std::max(0.0, infected_distance_square_sum)
                / double(infected.size())
            );
        }
        snapshot.max_ever_distance = max_ever_distance;
        const int shell_limit = std::min(
            int(shell_q_count.size()),
            std::max(cfg.information_radius + 3, 9)
        );
        snapshot.shell_q_count.assign(
            shell_q_count.begin(), shell_q_count.begin() + shell_limit
        );
        detail_bins[bin] = snapshot;
    }

    std::string serialize_population_interface_trace() const {
        std::vector<int> bins;
        bins.reserve(detail_bins.size());
        for (const auto &[index, snapshot] : detail_bins) {
            (void)snapshot;
            bins.push_back(index);
        }
        std::sort(bins.begin(), bins.end());
        std::ostringstream stream;
        for (size_t position = 0; position < bins.size(); ++position) {
            if (position) stream << ';';
            const int index = bins[position];
            const auto &snapshot = detail_bins.at(index);
            stream << index;
            for (int state_index = 0; state_index < 5; ++state_index) {
                stream << ':' << snapshot.state_count[state_index];
            }
            stream << ':' << snapshot.si_n_edges
                   << ':' << snapshot.si_q_edges
                   << ':' << snapshot.q_interface_nodes
                   << ':' << snapshot.q_far_nodes
                   << ':' << snapshot.infected_radius_mean
                   << ':' << snapshot.infected_radius_rms
                   << ':' << snapshot.max_ever_distance;
        }
        return stream.str();
    }

    std::string serialize_shell_trace() const {
        std::vector<int> bins;
        bins.reserve(detail_bins.size());
        for (const auto &[index, snapshot] : detail_bins) {
            (void)snapshot;
            bins.push_back(index);
        }
        std::sort(bins.begin(), bins.end());
        std::ostringstream stream;
        for (size_t position = 0; position < bins.size(); ++position) {
            if (position) stream << ';';
            const int index = bins[position];
            const auto &snapshot = detail_bins.at(index);
            stream << index << ':';
            for (size_t shell = 0; shell < snapshot.shell_q_count.size(); ++shell) {
                if (shell) stream << ',';
                stream << snapshot.shell_q_count[shell]
                       << '/' << shell_total_count[shell];
            }
        }
        return stream.str();
    }

    void initial_geometry(int seed, Outcome &out) const {
        if (!cfg.instrumentation) return;
        const int n = graph.n;
        std::vector<char> is_q_node(n, 0);
        int q_total = 0;
        for (int u = 0; u < n; ++u) {
            is_q_node[u] = is_q(state[u]);
            q_total += is_q_node[u];
        }
        out.initial_q = q_total;
        out.min_q = q_total;
        if (q_total == 0) {
            out.initial_q_components = 0;
            out.initial_q_largest_component = 0;
            out.initial_q_largest_fraction = 0.0;
            out.initial_q_component_second_moment = 0.0;
            out.initial_qn_interface = 0;
            out.initial_qn_interface_fraction = 0.0;
            out.seed_distance_to_q = -1;
            out.seed_q_coverage_r1 = 0.0;
            out.seed_q_coverage_r2 = 0.0;
            out.seed_q_coverage_r4 = 0.0;
            return;
        }

        std::vector<char> visited(n, 0);
        std::vector<int> queue;
        queue.reserve(n);
        int components = 0;
        int largest = 0;
        double second_moment_sum = 0.0;
        for (int root = 0; root < n; ++root) {
            if (!is_q_node[root] || visited[root]) continue;
            components++;
            queue.clear();
            queue.push_back(root);
            visited[root] = 1;
            size_t head = 0;
            while (head < queue.size()) {
                const int u = queue[head++];
                for (int v : graph.neighbors[u]) {
                    if (is_q_node[v] && !visited[v]) {
                        visited[v] = 1;
                        queue.push_back(v);
                    }
                }
            }
            const int size = int(queue.size());
            largest = std::max(largest, size);
            second_moment_sum += double(size) * size;
        }
        out.initial_q_components = components;
        out.initial_q_largest_component = largest;
        out.initial_q_largest_fraction = double(largest) / q_total;
        out.initial_q_component_second_moment = second_moment_sum / q_total;

        int interface_count = 0;
        for (const auto &[u, v] : graph.edges) {
            interface_count += is_q_node[u] != is_q_node[v];
        }
        out.initial_qn_interface = interface_count;
        out.initial_qn_interface_fraction =
            double(interface_count) / double(graph.edges.size());

        std::vector<int> distance(n, -1);
        queue.clear();
        queue.push_back(seed);
        distance[seed] = 0;
        size_t head = 0;
        std::array<int, 5> ball_total{};
        std::array<int, 5> ball_q{};
        while (head < queue.size()) {
            const int u = queue[head++];
            const int d = distance[u];
            if (d > 4) continue;
            if (d > 0) {
                for (int radius : {1, 2, 4}) {
                    if (d <= radius) {
                        ball_total[radius]++;
                        ball_q[radius] += is_q_node[u];
                    }
                }
            }
            if (is_q_node[u] && out.seed_distance_to_q < 0) {
                out.seed_distance_to_q = d;
            }
            if (d == 4) continue;
            for (int v : graph.neighbors[u]) {
                if (distance[v] < 0) {
                    distance[v] = d + 1;
                    queue.push_back(v);
                }
            }
        }
        auto coverage = [&](int radius) {
            return ball_total[radius] > 0
                ? double(ball_q[radius]) / ball_total[radius]
                : std::numeric_limits<double>::quiet_NaN();
        };
        out.seed_q_coverage_r1 = coverage(1);
        out.seed_q_coverage_r2 = coverage(2);
        out.seed_q_coverage_r4 = coverage(4);
    }

    static void synchronization_summary(
        const std::unordered_map<int, int> &bins,
        int total_events,
        double bin_width,
        int population,
        double infection_peak_time,
        int &max_bin,
        double &max_avalanche,
        double &concentration,
        double &effective_bins,
        double &peak_time,
        double &lag
    ) {
        if (total_events <= 0) return;
        int peak_index = -1;
        int peak_count = -1;
        double square_sum = 0.0;
        for (const auto &[index, count] : bins) {
            square_sum += double(count) * count;
            if (count > peak_count || (count == peak_count && index < peak_index)) {
                peak_count = count;
                peak_index = index;
            }
        }
        max_bin = peak_count;
        max_avalanche = double(peak_count) / population;
        concentration = square_sum / (double(total_events) * total_events);
        effective_bins = 1.0 / concentration;
        peak_time = (peak_index + 0.5) * bin_width;
        lag = peak_time - infection_peak_time;
    }

    static std::string serialize_histogram(
        const std::unordered_map<int, int> &bins
    ) {
        std::vector<std::pair<int, int>> ordered(bins.begin(), bins.end());
        std::sort(ordered.begin(), ordered.end());
        std::ostringstream stream;
        for (size_t index = 0; index < ordered.size(); ++index) {
            if (index) stream << ';';
            stream << ordered[index].first << ':' << ordered[index].second;
        }
        return stream.str();
    }

    std::string serialize_iq_trace() const {
        std::vector<int> bins;
        bins.reserve(i_bin_max.size());
        for (const auto &[index, value] : i_bin_max) {
            (void)value;
            bins.push_back(index);
        }
        std::sort(bins.begin(), bins.end());
        std::ostringstream stream;
        for (size_t position = 0; position < bins.size(); ++position) {
            if (position) stream << ';';
            const int index = bins[position];
            stream << index << ':' << i_bin_max.at(index) << ':' << q_bin_last.at(index);
        }
        return stream.str();
    }

    void finalize_instrumentation(
        Outcome &out,
        double start_i,
        double start_q,
        double last_i,
        double last_q,
        double twice_area
    ) const {
        out.max_ever_distance = max_ever_distance;
        out.rms_ever_distance = ever_distance_count > 0
            ? std::sqrt(ever_distance_square_sum / ever_distance_count)
            : 0.0;
        if (!cfg.instrumentation) return;
        out.trajectory_signed_area =
            0.5 * (twice_area + last_i * start_q - start_i * last_q);
        if (out.q_at_first_ic >= 0) {
            out.q_recovery_after_alert =
                double(std::max(0, out.max_q_after_alert - out.q_at_first_ic))
                / graph.n;
        }
        synchronization_summary(
            qn_bins, out.q_to_n_events, cfg.behavior_bin_width, graph.n,
            out.peak_time, out.qn_max_bin, out.qn_max_avalanche,
            out.qn_concentration, out.qn_effective_bins, out.qn_peak_time,
            out.qn_peak_lag_from_i_peak
        );
        synchronization_summary(
            nq_bins, out.n_to_q_events, cfg.behavior_bin_width, graph.n,
            out.peak_time, out.nq_max_bin, out.nq_max_avalanche,
            out.nq_concentration, out.nq_effective_bins, out.nq_peak_time,
            out.nq_peak_lag_from_i_peak
        );
        std::unordered_map<int, int> combined = nq_bins;
        for (const auto &[index, count] : qn_bins) combined[index] -= count;
        int peak_index = -1;
        int peak_value = 0;
        for (const auto &[index, value] : combined) {
            if (std::abs(value) > std::abs(peak_value)
                || (std::abs(value) == std::abs(peak_value)
                    && (peak_index < 0 || index < peak_index))) {
                peak_index = index;
                peak_value = value;
            }
        }
        if (peak_index >= 0) {
            out.net_flux_peak = double(peak_value) / graph.n;
            out.net_flux_peak_time = (peak_index + 0.5) * cfg.behavior_bin_width;
        }
        if (cfg.trace_mode) {
            out.qn_histogram = serialize_histogram(qn_bins);
            out.nq_histogram = serialize_histogram(nq_bins);
            out.iq_trace = serialize_iq_trace();
            out.population_interface_trace =
                serialize_population_interface_trace();
            out.shell_trace = serialize_shell_trace();
        }
    }

public:
    GraphProcess(
        const Config &config,
        const Graph &g,
        const GraphRuntime &runtime_state,
        const std::vector<std::vector<int>> *neighborhood,
        EdgePoolBackend backend = EdgePoolBackend::FenwickExact
    )
        : cfg(config), graph(g), runtime(runtime_state), edge_backend(backend), state(g.n, SN),
          adjacent_infected_count(g.n, 0),
          signal_neighborhood(neighborhood), infected(g.n),
          ever_infected(g.n, 0), q_rng(config.q_init_seed),
          i_rng(config.i_init_seed), rng(config.dynamic_seed) {
        const int directed = int(runtime.directed_u.size());
        infection_edge_category.assign(directed, -1);
        behavior_edge_category.assign(directed, -1);
        if (edge_backend == EdgePoolBackend::UniformRegular && !runtime.regular)
            throw std::runtime_error("uniform_regular requires a regular graph");
        for (int category = 0; category < 8; ++category) {
            const double uniform_weight = (category < 4 && cfg.contact_mode == "per_edge")
                ? 1.0 : (runtime.regular ? 1.0 / double(runtime.regular_degree) : 0.0);
            pools[category] = EdgePool(directed, edge_backend, uniform_weight);
        }
        if (uses_node_behavior_pools()) {
            for (auto &pool : behavior_node_pools) pool = Fenwick(graph.n);
            behavior_node_category.assign(graph.n, -1);
            behavior_node_weight.assign(graph.n, 0.0);
        }
        interface_affected_scratch.reserve(runtime.regular ? runtime.regular_degree + 1 : 16);
        local_affected_scratch.reserve(runtime.regular ? 2 * runtime.regular_degree + 8 : 32);
    }

    Outcome run() {
        const int initial_i = std::clamp(
            cfg.i0 > 0.0 ? int(std::llround(cfg.i0 * graph.n)) : 1, 1, graph.n
        );
        std::vector<int> q_order(graph.n);
        std::vector<int> i_order(graph.n);
        std::iota(q_order.begin(), q_order.end(), 0);
        std::iota(i_order.begin(), i_order.end(), 0);
        q_rng.shuffle(q_order);
        i_rng.shuffle(i_order);
        const int seed = i_order.front();
        prepare_seed_distances(seed);
        if (cfg.q0 > 0.0 && cfg.q_placement == "none") {
            const int q_total = std::clamp(int(std::llround(cfg.q0 * graph.n)), 0, graph.n);
            for (int i = 0; i < q_total; ++i) state[q_order[i]] = SQ;
        }
        const int q_seed = select_q_seed(seed);
        if (q_seed >= 0) state[q_seed] = SQ;
        for (int j = 0; j < initial_i; ++j) {
            const int u = i_order[j];
            state[u] = is_q(state[u]) ? IQ : IN;
        }
        for (int u = 0; u < graph.n; ++u) {
            if (is_q(state[u]) && distance_from_seed[u] >= 0) {
                shell_q_count[distance_from_seed[u]]++;
            }
        }
        for (int s : state) count[s]++;
        infected_distance_sum = 0.0;
        infected_distance_square_sum = 0.0;
        for (int u = 0; u < graph.n; ++u) {
            if (is_infected(state[u])) { infected.add(u); mark_ever(u); }
        }
        for (int e = 0; e < int(runtime.directed_u.size()); ++e) refresh(e);
        initialize_information_and_interface();

        Outcome out;
        out.initial_i = initial_i;
        out.initial_iq = int(count[IQ]);
        out.ever_infected_count = ever_infected_count;
        out.signal_neighborhood_size =
            signal_neighborhood == nullptr || signal_neighborhood->empty()
                ? 0 : int((*signal_neighborhood)[seed].size());
        out.seed_vertex = seed;
        out.seed_degree = graph.degree[seed];
        out.q_seed_vertex = q_seed;
        out.q_seed_degree = q_seed >= 0 ? graph.degree[q_seed] : 0;
        out.q_seed_distance = q_seed < 0 ? -1 : (cfg.q_placement == "adjacent" ? 1 : -2);
        out.peak_infected = int(count[IN] + count[IQ]);
        out.max_q = int(count[SQ] + count[IQ]);
        out.q_mesoscopic = out.max_q >= int(std::ceil(std::sqrt(double(graph.n))));
        out.q_macroscopic = double(out.max_q) / graph.n >= 0.05;
        out.activated_ic = double(infected.size()) / graph.n >= cfg.i_c;
        initial_geometry(seed, out);
        if (!cfg.instrumentation) {
            out.initial_q = int(count[SQ] + count[IQ]);
            out.min_q = out.initial_q;
        }
        if (out.activated_ic && cfg.instrumentation) {
            out.first_ic_time = 0.0;
            out.q_at_first_ic = out.initial_q;
            out.initial_q_lost_before_ic = 0.0;
            out.max_q_after_alert = out.initial_q;
        }
        if (cfg.instrumentation && cfg.trace_mode) {
            i_bin_max[0] = initial_i;
            q_bin_last[0] = out.initial_q;
            capture_detail_bin(0);
        }
        initialize_trace_v2();
        double time = 0.0;
        const double start_i = double(infected.size()) / graph.n;
        const double start_q = double(out.initial_q) / graph.n;
        double previous_i_fraction = start_i;
        double previous_q_fraction = start_q;
        double twice_area = 0.0;
        const double beta_n = cfg.beta_n;
        const double beta_q = cfg.beta_q;
        const double beta_a = cfg.beta_a;
        while (infected.size() > 0
               && out.accepted_event_count < cfg.event_limit
               && time < cfg.tmax) {
            const int current_i = int(infected.size());
            const bool node_information = uses_node_behavior_pools();
            const bool mixed_information = uses_mixed_thinning();
            std::pair<double, double> p{0.0, 0.0};
            if (!node_information && !mixed_information
                && cfg.behavior_mode == "adaptive") {
                p = global_adoption_cached(current_i);
            }
            std::array<double, 9> rates{
                beta_n * pools[0].total(),
                beta_a * pools[1].total(),
                beta_a * pools[2].total(),
                beta_q * pools[3].total(),
                node_information ? behavior_node_pools[0].total()
                    : cfg.tau * (mixed_information ? 1.0 : p.first)
                        * pools[4].total(),
                node_information ? behavior_node_pools[1].total()
                    : cfg.tau * (mixed_information ? 1.0 : p.second)
                        * pools[5].total(),
                node_information ? behavior_node_pools[2].total()
                    : cfg.tau * (mixed_information ? 1.0 : p.first)
                        * pools[6].total(),
                node_information ? behavior_node_pools[3].total()
                    : cfg.tau * (mixed_information ? 1.0 : p.second)
                        * pools[7].total(),
                double(current_i)
            };
            const double total = std::accumulate(rates.begin(), rates.end(), 0.0);
            if (!(total > 0.0) || !std::isfinite(total)) break;
            const double previous_time = time;
            const double next_time = time
                + -std::log(std::max(1e-300, rng.uniform01())) / total;
            fill_trace_v2(
                std::min(next_time, cfg.tmax), current_i,
                int(count[SQ] + count[IQ])
            );
            time = std::min(next_time, cfg.tmax);
            if (cfg.instrumentation
                && double(current_i) / graph.n >= cfg.i_c) {
                out.time_above_ic += std::min(time, cfg.tmax) - previous_time;
            }
            const double interval = std::min(time, cfg.tmax) - previous_time;
            out.q_person_time += double(count[SQ] + count[IQ]) / graph.n * interval;
            out.q_interface_person_time +=
                double(q_interface_nodes) / graph.n * interval;
            out.q_far_person_time +=
                double(count[SQ] + count[IQ] - q_interface_nodes)
                / graph.n * interval;
            out.si_n_contact_time +=
                double(raw_pool_counts[0] + raw_pool_counts[1])
                / graph.edges.size() * interval;
            out.si_q_contact_time +=
                double(raw_pool_counts[2] + raw_pool_counts[3])
                / graph.edges.size() * interval;
            out.n_infection_hazard_time +=
                (beta_n * pools[0].total() + beta_a * pools[1].total()) * interval;
            out.q_infection_hazard_time +=
                (beta_a * pools[2].total() + beta_q * pools[3].total()) * interval;
            if (next_time > cfg.tmax) break;
            const int previous_q = int(count[SQ] + count[IQ]);
            const bool previously_above_ic = double(current_i) / graph.n >= cfg.i_c;
            const int event = choose_event(rates, rng.uniform01() * total);
            double focal_local_signal = std::numeric_limits<double>::quiet_NaN();
            bool accepted_event = true;
            if (event < 8) {
                int focal = -1;
                int model = -1;
                if (event < 4 || !node_information) {
                    const double pool_total = pools[event].total();
                    const int edge = pools[event].sample(
                        std::min(std::nextafter(pool_total, 0.0), rng.uniform01() * pool_total)
                    );
                    focal = runtime.directed_u[edge];
                    model = runtime.directed_v[edge];
                } else {
                    const int local_category = event - 4;
                    const double pool_total = behavior_node_pools[local_category].total();
                    focal = behavior_node_pools[local_category].sample(
                        std::min(std::nextafter(pool_total, 0.0), rng.uniform01() * pool_total)
                    );
                }
                focal_local_signal = focal_local_prevalence(focal);
                if (mixed_information && event >= 4) {
                    const bool n_to_q = event == 4 || event == 6;
                    // Pairwise imitation: N->Q depends on the focal N payoff;
                    // Q->N depends on the copied N model's own information/payoff.
                    const int signal_owner = n_to_q ? focal : model;
                    const auto mixed_p =
                        adoption_from_prevalence(cfg, information_signal(signal_owner));
                    const double acceptance =
                        n_to_q ? mixed_p.first : mixed_p.second;
                    accepted_event = rng.uniform01() < acceptance;
                    if (!accepted_event) out.null_behavior_proposals++;
                }
                if (!accepted_event) {
                    out.events++;
                    if (cfg.instrumentation && cfg.trace_mode) {
                        const int bin = int(std::floor(time / cfg.behavior_bin_width));
                        capture_detail_bin(bin);
                    }
                    continue;
                }
                if (event == 0 || event == 1) out.infections_into_n++;
                if (event == 2 || event == 3) out.infections_into_q++;
                if (event == 5 || event == 7) {
                    out.qn_local_prevalence_sum += focal_local_signal;
                    out.qn_near_infected_events += focal_local_signal > 0.0;
                } else if (event == 4 || event == 6) {
                    out.nq_local_prevalence_sum += focal_local_signal;
                    out.nq_near_infected_events += focal_local_signal > 0.0;
                }
                switch (event) {
                    case 0: case 1: change_state(focal, IN); break;
                    case 2: case 3: change_state(focal, IQ); break;
                    case 4: change_state(focal, SQ); break;
                    case 5: change_state(focal, SN); break;
                    case 6: change_state(focal, IQ); break;
                    case 7: change_state(focal, IN); break;
                }
            } else {
                const int focal = infected.sample(rng);
                change_state(focal, RR);
            }
            out.events++;
            out.accepted_event_count++;
            const int now_i = int(infected.size());
            const int now_q = int(count[SQ] + count[IQ]);
            if (cfg.instrumentation) {
                const int bin = int(std::floor(time / cfg.behavior_bin_width));
                if (event == 5 || event == 7) {
                    out.q_to_n_events++;
                    qn_bins[bin]++;
                } else if (event == 4 || event == 6) {
                    out.n_to_q_events++;
                    nq_bins[bin]++;
                }
                if (cfg.trace_mode) {
                    auto found = i_bin_max.find(bin);
                    if (found == i_bin_max.end()) {
                        i_bin_max[bin] = now_i;
                    } else {
                        found->second = std::max(found->second, now_i);
                    }
                    q_bin_last[bin] = now_q;
                    capture_detail_bin(bin);
                }
                out.min_q = std::min(out.min_q, now_q);
                const bool now_above_ic = double(now_i) / graph.n >= cfg.i_c;
                if (now_above_ic != previously_above_ic) out.ic_crossings++;
                if (!previously_above_ic && now_above_ic
                    && !std::isfinite(out.first_ic_time)) {
                    out.first_ic_time = time;
                    out.q_at_first_ic = now_q;
                    out.initial_q_lost_before_ic = out.initial_q > 0
                        ? double(std::max(0, out.initial_q - now_q)) / out.initial_q
                        : std::numeric_limits<double>::quiet_NaN();
                    out.max_q_after_alert = now_q;
                }
                if (std::isfinite(out.first_ic_time)) {
                    out.max_q_after_alert = std::max(out.max_q_after_alert, now_q);
                }
                const double now_i_fraction = double(now_i) / graph.n;
                const double now_q_fraction = double(now_q) / graph.n;
                twice_area += previous_i_fraction * now_q_fraction
                    - now_i_fraction * previous_q_fraction;
                previous_i_fraction = now_i_fraction;
                previous_q_fraction = now_q_fraction;
            }
            if (now_i > out.peak_infected) {
                out.peak_infected = now_i;
                out.peak_time = time;
            }
            out.max_q = std::max(out.max_q, int(count[SQ] + count[IQ]));
            if (previous_q > 0 && count[SQ] + count[IQ] == 0) out.q_extinction_time = time;
            out.q_mesoscopic = out.q_mesoscopic
                || out.max_q >= int(std::ceil(std::sqrt(double(graph.n))));
            out.q_macroscopic = out.q_macroscopic
                || double(out.max_q) / graph.n >= 0.05;
            out.activated_ic = out.activated_ic || double(now_i) / graph.n >= cfg.i_c;
        }
        out.censored = infected.size() > 0;
        out.ever_infected_count = ever_infected_count;
        out.outbreak_size = ever_infected_count;
        out.duration = time;
        out.stop_reason = out.censored
            ? (time >= cfg.tmax ? "tmax" : "event_limit")
            : "absorbed_I0";
        out.final_q = int(count[SQ] + count[IQ]);
        if (out.max_q > 0 && out.final_q == 0 && !std::isfinite(out.q_extinction_time)) {
            out.q_extinction_time = time;
        }
        out.infected_degree_mean = out.outbreak_size > 0
            ? double(infected_degree_sum) / double(out.outbreak_size)
            : std::numeric_limits<double>::quiet_NaN();
        out.hubs_infected = hubs_infected;
        const double total_hazard_time =
            out.n_infection_hazard_time + out.q_infection_hazard_time;
        if (total_hazard_time > 0.0) {
            out.q_hazard_fraction = out.q_infection_hazard_time / total_hazard_time;
        }
        if (out.q_to_n_events > 0) {
            out.qn_mean_local_prevalence =
                out.qn_local_prevalence_sum / out.q_to_n_events;
        }
        if (out.n_to_q_events > 0) {
            out.nq_mean_local_prevalence =
                out.nq_local_prevalence_sum / out.n_to_q_events;
        }
        finalize_instrumentation(
            out, start_i, start_q, previous_i_fraction, previous_q_fraction, twice_area
        );
        out.trace_v2 = serialize_trace_v2(
            time, int(infected.size()), int(count[SQ] + count[IQ])
        );
        return out;
    }
};

// Strict FPM_RUN_SPEC_V6 parser. Fail closed on missing/extra fields, invalid
// domains, or seed collisions.  This keeps configuration drift from becoming a
// silent change in scientific semantics.
static std::vector<Config> read_design(const std::string &path, bool allow_legacy_identity, bool allow_legacy_beta_hierarchy) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("could not open design");
    std::string header_line;
    std::getline(input, header_line);
    const auto header = split(header_line, '\t');
    const std::vector<std::string> expected_header{
        "schema_version", "analysis_role", "run_fingerprint", "run_index", "run_id",
        "stage", "topology", "graph_path", "graph_sha256", "graph_certificate_path",
        "graph_certificate_sha256", "boundary_condition", "population", "n0",
        "n_occ", "n_gc", "graph_seed", "q_init_seed", "i_init_seed", "dynamic_seed",
        "contact_mode", "beta_n", "beta_a", "beta_q", "beta_mix_a", "beta_mixing_rule", "profile", "arm", "protocol", "q0", "i0",
        "tau", "i_c", "fermi_k", "omega", "delta", "seed_policy",
        "behavior_mode", "q_placement", "information_mode", "information_radius",
        "information_alpha", "signal_seed", "instrumentation_mode", "trace_mode",
        "behavior_bin_width", "event_limit", "tmax"
    };
    auto expected_header_extended = expected_header;
    expected_header_extended.push_back("fingerprint_version");
    auto v5_header = expected_header;
    for (auto &name : v5_header) if (name == "fermi_k") name = "epsilon";
    auto v5_header_extended = v5_header;
    v5_header_extended.push_back("fingerprint_version");
    auto v4_header = v5_header;
    for (auto &name : v4_header) if (name == "tau") name = "chi";
    auto v4_header_extended = v4_header;
    v4_header_extended.push_back("fingerprint_version");
    const bool current_extended = header == expected_header_extended;
    const bool legacy_v5_extended = header == v5_header_extended;
    const bool legacy_v4_extended = header == v4_header_extended;
    const bool legacy_v5 = header == v5_header || legacy_v5_extended;
    const bool legacy_v4 = header == v4_header || legacy_v4_extended;
    const bool extended_identity = current_extended || legacy_v5_extended || legacy_v4_extended;
    if (header != expected_header && !current_extended && !legacy_v5 && !legacy_v4)
        throw std::runtime_error("design schema/header mismatch");
    if (!extended_identity && !allow_legacy_identity)
        throw std::runtime_error("opaque run identity requires --allow-legacy-identity");
    std::unordered_map<std::string, int> column;
    for (int i = 0; i < int(header.size()); ++i) {
        if (!column.emplace(header[i], i).second) {
            throw std::runtime_error("duplicate design column: " + header[i]);
        }
    }
    auto value = [&](const std::vector<std::string> &fields, const std::string &name)
        -> const std::string & { return fields.at(column.at(name)); };
    auto optional = [&](const std::vector<std::string> &fields, const std::string &name,
                        const std::string &fallback) -> std::string {
        const auto found = column.find(name);
        return found == column.end() ? fallback : fields.at(found->second);
    };
    std::vector<Config> configs;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto f = split(line, '\t');
        if (f.size() != header.size()) throw std::runtime_error("design row width mismatch");
        Config c;
        c.schema_version = value(f, "schema_version");
        c.analysis_role = value(f, "analysis_role");
        c.run_fingerprint = value(f, "run_fingerprint");
        c.fingerprint_version = optional(f, "fingerprint_version", "");
        if (!is_lower_hex_sha256(c.run_fingerprint)
            || !is_lower_hex_sha256(value(f, "graph_sha256"))
            || !is_lower_hex_sha256(value(f, "graph_certificate_sha256")))
            throw std::runtime_error("run, graph and certificate identities must be lowercase hexadecimal SHA-256");
        if (extended_identity) {
            const std::string expected_version = legacy_v4_extended ? "FPM_TSV_ID_V1"
                : legacy_v5_extended ? "FPM_TSV_ID_V2" : "FPM_TSV_ID_V3";
            if (c.fingerprint_version != expected_version)
                throw std::runtime_error("unsupported fingerprint_version");
            const auto expected_fp = fpm::sha256_hex(canonical_semantic_identity_json(f, column));
            if (c.run_fingerprint != expected_fp)
                throw std::runtime_error("semantic run_fingerprint mismatch");
        }
        c.run_index = std::stoi(value(f, "run_index"));
        c.run_id = value(f, "run_id");
        c.stage = value(f, "stage");
        c.topology = value(f, "topology");
        c.graph_path = value(f, "graph_path");
        c.graph_sha256 = value(f, "graph_sha256");
        c.graph_certificate_path = value(f, "graph_certificate_path");
        c.graph_certificate_sha256 = value(f, "graph_certificate_sha256");
        c.boundary_condition = value(f, "boundary_condition");
        c.population = std::stoi(value(f, "population"));
        c.n0 = std::stoi(optional(f, "n0", std::to_string(c.population)));
        c.n_occ = std::stoi(optional(f, "n_occ", std::to_string(c.population)));
        c.n_gc = std::stoi(optional(f, "n_gc", std::to_string(c.population)));
        c.graph_seed = std::stoull(value(f, "graph_seed"));
        c.q_init_seed = std::stoull(value(f, "q_init_seed"));
        c.i_init_seed = std::stoull(value(f, "i_init_seed"));
        c.dynamic_seed = std::stoull(value(f, "dynamic_seed"));
        c.contact_mode = value(f, "contact_mode");
        c.beta_n = std::stod(value(f, "beta_n"));
        c.profile = value(f, "profile");
        c.arm = optional(f, "arm", c.profile);
        c.protocol = optional(f, "protocol", "LEGACY_IC_FIXED");
        c.q0 = std::stod(value(f, "q0"));
        c.i0 = std::stod(optional(f, "i0", "0.0"));
        c.beta_q = std::stod(value(f, "beta_q"));
        c.beta_mix_a = std::stod(value(f, "beta_mix_a"));
        c.beta_mixing_rule = value(f, "beta_mixing_rule");
        const auto resolved = fpm::resolve_betas(c.beta_n, c.beta_q, c.beta_mix_a, fpm::parse_beta_mixing_rule(c.beta_mixing_rule), allow_legacy_beta_hierarchy);
        c.beta_a = resolved.beta_a;
        c.claim_grade_corrected_model = resolved.claim_grade_corrected_model;
        const double declared_beta_a = std::stod(value(f, "beta_a"));
        if (!std::isfinite(declared_beta_a) || std::abs(declared_beta_a-c.beta_a) > 1e-12*std::max(1.0,std::abs(c.beta_a)))
            throw std::runtime_error("declared beta_a does not match centralized resolver");
        c.tau = std::stod(value(f, legacy_v4 ? "chi" : "tau"));
        c.i_c = std::stod(value(f, "i_c"));
        c.fermi_k = std::stod(value(f, (legacy_v4 || legacy_v5) ? "epsilon" : "fermi_k"));
        c.omega = std::stod(optional(f, "omega", "1.0"));
        c.delta = std::stod(optional(
            f, "delta", std::to_string(c.omega / (c.i_c * c.beta_n))
        ));
        c.seed_policy = value(f, "seed_policy");
        c.behavior_mode = optional(
            f, "behavior_mode", c.profile == "EGT" ? "adaptive" : "off"
        );
        c.q_placement = optional(f, "q_placement", "none");
        c.information_mode = optional(f, "information_mode", "global");
        c.information_radius = std::stoi(
            optional(f, "information_radius", "1")
        );
        c.information_alpha = std::stod(
            optional(f, "information_alpha", "0")
        );
        c.signal_seed = std::stoull(optional(
            f, "signal_seed", std::to_string(c.dynamic_seed ^ 0x9e3779b97f4a7c15ULL)
        ));
        c.instrumentation = optional(f, "instrumentation_mode", "on") != "off";
        c.trace_mode = optional(f, "trace_mode", "off") == "on";
        c.behavior_bin_width = std::stod(
            optional(f, "behavior_bin_width", "0.25")
        );
        c.event_limit = std::stoll(value(f, "event_limit"));
        c.tmax = std::stod(value(f, "tmax"));
        const auto finite = [](double value) { return std::isfinite(value); };
        const bool spatial_periodic = fpm::is_spatial_periodic_topology(c.topology);
        const std::string zero_sha256(64, '0');
        const bool legacy_certificate_bypass =
            spatial_periodic && allow_legacy_identity && !extended_identity;
        c.claim_grade_materiality = !legacy_certificate_bypass;
        if ((legacy_v4 ? c.schema_version != "FPM_RUN_SPEC_V4"
                : legacy_v5 ? c.schema_version != "FPM_RUN_SPEC_V5"
                : c.schema_version != "FPM_RUN_SPEC_V6")
            || c.population < 2 || !finite(c.beta_n) || c.beta_n <= 0
            || !finite(c.q0) || c.q0 < 0 || c.q0 > 1
            || !finite(c.i_c) || c.i_c <= 0
            || !finite(c.fermi_k) || c.fermi_k <= 0
            || !finite(c.omega) || c.omega <= 0
            || !finite(c.delta) || c.delta <= 0
            || !finite(c.i0) || c.i0 < 0 || c.i0 > 1
            || !finite(c.beta_q) || c.beta_q <= 0
            || !finite(c.beta_a) || !finite(c.beta_mix_a)
            || !finite(c.tau) || c.tau < 0
            || c.event_limit <= 0 || c.tmax <= 0 || c.n0 < c.n_occ
            || c.n_occ < c.n_gc || c.n_gc != c.population
            || !finite(c.behavior_bin_width) || c.behavior_bin_width <= 0
            || !finite(c.tmax) || c.information_radius < 1
            || c.information_radius > 16 || !finite(c.information_alpha)
            || c.information_alpha < 0.0
            || c.information_alpha > 1.0
            || c.q_init_seed == c.i_init_seed
            || c.q_init_seed == c.dynamic_seed
            || c.i_init_seed == c.dynamic_seed
            || c.signal_seed == c.q_init_seed
            || c.signal_seed == c.i_init_seed
            || c.signal_seed == c.dynamic_seed
            || (spatial_periodic && !legacy_certificate_bypass
                && (c.graph_certificate_path != c.graph_path + ".periodic.json"
                    || c.graph_certificate_sha256 == zero_sha256
                    || c.boundary_condition != "periodic_torus"))
            || ((!spatial_periodic || legacy_certificate_bypass)
                && (c.graph_certificate_path != "-"
                    || c.graph_certificate_sha256 != zero_sha256
                    || c.boundary_condition != (legacy_certificate_bypass
                        ? "legacy_unverified"
                        : (c.topology == "well_mixed" ? "well_mixed" : "nonspatial_network"))))) {
            throw std::runtime_error("invalid config");
        }
        static const std::unordered_set<std::string> roles{
            "primary", "trace_diagnostic", "benchmark", "calibration", "smoke_test"
        };
        if (!roles.count(c.analysis_role)) throw std::runtime_error("invalid analysis role");
        if (c.behavior_mode != "frozen" && c.behavior_mode != "adaptive") {
            throw std::runtime_error("invalid behavior mode");
        }
        if (c.contact_mode != "normalized") throw std::runtime_error("unsupported contact mode");
        if (c.seed_policy != "uniform") throw std::runtime_error("unsupported seed policy");
        if (c.q_placement != "none") throw std::runtime_error("unsupported Q placement");
        if (c.information_mode != "global" && c.information_mode != "local" && c.information_mode != "local_shuffled" && c.information_mode != "mixed") throw std::runtime_error("unsupported information mode");
        if (c.trace_mode != (c.analysis_role == "trace_diagnostic")) {
            throw std::runtime_error("trace role/mode mismatch");
        }
        configs.push_back(c);
    }
    return configs;
}

// Generic CLI entry point.  Topology, mission/stage and parameter regime are all
// data in DESIGN.tsv; no topology- or mission-specific executables are required.
int main(int argc, char **argv) {
    if (argc < 3 || argc > 6) {
        std::cerr << "usage: fpm_network_sim DESIGN.tsv SUMMARY.tsv [fenwick_exact|uniform_regular] [--allow-legacy-identity] [--allow-legacy-beta-hierarchy]\n";
        return 2;
    }
    EdgePoolBackend backend = EdgePoolBackend::FenwickExact;
    bool backend_seen=false, allow_legacy_identity=false, allow_legacy_beta_hierarchy=false;
    for (int i=3;i<argc;++i) {
        const std::string arg=argv[i];
        if (arg=="fenwick_exact" || arg=="uniform_regular") {
            if (backend_seen) { std::cerr << "duplicate edge-pool backend\n"; return 2; }
            backend_seen=true;
            backend = arg=="uniform_regular" ? EdgePoolBackend::UniformRegular : EdgePoolBackend::FenwickExact;
        } else if (arg=="--allow-legacy-identity") {
            allow_legacy_identity=true;
        } else if (arg=="--allow-legacy-beta-hierarchy") {
            allow_legacy_beta_hierarchy=true;
        } else {
            std::cerr << "unknown argument: " << arg << '\n'; return 2;
        }
    }
    try {
        const auto configs = read_design(argv[1], allow_legacy_identity, allow_legacy_beta_hierarchy);
        std::ofstream output(argv[2]);
        if (!output) throw std::runtime_error("could not open summary output");
        output << "schema_version\tanalysis_role\trun_fingerprint"
               << "\trun_index\trun_id\tstage\ttopology\tpopulation"
               << "\tgraph_sha256\tgraph_certificate_sha256\tboundary_condition"
               << "\tclaim_grade_materiality"
               << "\tgraph_seed\tq_init_seed\ti_init_seed"
               << "\tdynamic_seed\tcontact_mode\tprotocol\tomega\tdelta"
               << "\tbeta_n\tbeta_a\tbeta_q\tbeta_mix_a\tbeta_mixing_rule\tclaim_grade_corrected_model\ttau\ti_c"
               << "\tfermi_k\tprofile\tseed_policy"
               << "\tevents\taccepted_event_count\tcensored\tstop_reason"
               << "\toutbreak_size\tever_infected_count\tduration\tpeak_infected"
               << "\tpeak_time\tmax_q\tfinal_q\tq_extinction_time\tq_person_time"
               << "\tq_interface_person_time\tq_far_person_time"
               << "\tsi_q_contact_time\tsi_n_contact_time"
               << "\tq_mesoscopic\tq_macroscopic\tactivated_ic\tseed_vertex"
               << "\tseed_degree\tq_seed_vertex\tq_seed_degree\tq_seed_distance"
               << "\tinfected_degree_mean\thubs_infected\tn0\tn_occ\tn_gc"
               << "\tbehavior_mode\tq_placement\tarm\tq0\ti0\tinformation_mode"
               << "\tinformation_radius\tinformation_alpha"
               << "\tsignal_seed\tinstrumentation_mode"
               << "\tbehavior_bin_width\tinitial_q\tinitial_i\tinitial_iq"
               << "\tmin_q\tfirst_ic_time"
               << "\tq_at_first_ic\tinitial_q_lost_before_ic\ttime_above_ic"
               << "\tic_crossings\tq_to_n_events\tn_to_q_events"
               << "\tinfections_into_n\tinfections_into_q"
               << "\tn_infection_hazard_time\tq_infection_hazard_time"
               << "\tq_hazard_fraction\tqn_local_prevalence_sum"
               << "\tnq_local_prevalence_sum\tqn_near_infected_events"
               << "\tnq_near_infected_events\tqn_mean_local_prevalence"
               << "\tnq_mean_local_prevalence"
               << "\tmax_q_after_alert\tq_recovery_after_alert"
               << "\ttrajectory_signed_area\tqn_max_bin\tnq_max_bin"
               << "\tqn_max_avalanche\tnq_max_avalanche"
               << "\tqn_concentration\tnq_concentration"
               << "\tqn_effective_bins\tnq_effective_bins"
               << "\tqn_peak_time\tnq_peak_time"
               << "\tqn_peak_lag_from_i_peak\tnq_peak_lag_from_i_peak"
               << "\tnet_flux_peak\tnet_flux_peak_time"
               << "\tinitial_q_components\tinitial_q_largest_component"
               << "\tinitial_q_largest_fraction"
               << "\tinitial_q_component_second_moment"
               << "\tinitial_qn_interface\tinitial_qn_interface_fraction"
               << "\tseed_distance_to_q\tseed_q_coverage_r1"
               << "\tseed_q_coverage_r2\tseed_q_coverage_r4"
               << "\tnull_behavior_proposals\tmax_ever_distance"
               << "\tsignal_neighborhood_size\trms_ever_distance"
               << "\ttrace_mode\tqn_histogram\tnq_histogram\tiq_trace"
               << "\tpopulation_interface_trace\tshell_trace\ttrace_v2\n";
        std::string cached_path;
        std::string cached_graph_sha256;
        std::string cached_graph_certificate_sha256;
        std::string cached_boundary_condition;
        std::string cached_topology;
        Graph graph;
        GraphRuntime graph_runtime;
        std::unordered_map<int, std::vector<std::vector<int>>> neighborhood_cache;
        for (const auto &cfg : configs) {
            Outcome out;
            if (cfg.topology == "well_mixed") {
                out = run_well_mixed(cfg);
            } else {
                if (cfg.graph_path != cached_path) {
                    cached_graph_sha256 = fpm::sha256_file(cfg.graph_path);
                    if (cached_graph_sha256 != cfg.graph_sha256)
                        throw std::runtime_error("graph SHA-256 mismatch");
                    if (fpm::is_spatial_periodic_topology(cfg.topology)
                        && cfg.claim_grade_materiality) {
                        const auto certificate = fpm::validate_periodic_certificate(
                            cfg.graph_certificate_path,
                            cfg.graph_certificate_sha256,
                            cfg.topology,
                            cached_graph_sha256
                        );
                        cached_graph_certificate_sha256 = certificate.sha256;
                    } else {
                        cached_graph_certificate_sha256 = cfg.graph_certificate_sha256;
                    }
                    cached_boundary_condition = cfg.boundary_condition;
                    cached_topology = cfg.topology;
                    graph = load_graph(cfg.graph_path, cfg.topology);
                    graph_runtime = build_graph_runtime(graph);
                    cached_path = cfg.graph_path;
                    neighborhood_cache.clear();
                } else if (cached_graph_sha256 != cfg.graph_sha256
                           || cached_graph_certificate_sha256 != cfg.graph_certificate_sha256
                           || cached_boundary_condition != cfg.boundary_condition
                           || cached_topology != cfg.topology) {
                    throw std::runtime_error("cached graph materiality mismatch");
                }
                if (graph.n != cfg.population) throw std::runtime_error("graph/config N mismatch");
                const std::vector<std::vector<int>> *neighborhood = nullptr;
                if (cfg.information_mode == "local"
                    || cfg.information_mode == "local_shuffled"
                    || cfg.information_mode == "mixed") {
                    auto found = neighborhood_cache.find(cfg.information_radius);
                    if (found == neighborhood_cache.end()) {
                        found = neighborhood_cache.emplace(
                            cfg.information_radius,
                            build_signal_neighborhood(graph, cfg.information_radius)
                        ).first;
                    }
                    neighborhood = &found->second;
                }
                GraphProcess process(cfg, graph, graph_runtime, neighborhood, backend);
                out = process.run();
            }
            output << cfg.schema_version << '\t' << cfg.analysis_role << '\t'
                   << cfg.run_fingerprint << '\t'
                   << cfg.run_index << '\t' << cfg.run_id << '\t' << cfg.stage << '\t'
                    << cfg.topology << '\t' << cfg.population
                    << '\t' << cfg.graph_sha256
                    << '\t' << cfg.graph_certificate_sha256
                    << '\t' << cfg.boundary_condition
                    << '\t' << (cfg.claim_grade_materiality ? 1 : 0)
                    << '\t' << cfg.graph_seed
                   << '\t' << cfg.q_init_seed << '\t' << cfg.i_init_seed
                   << '\t' << cfg.dynamic_seed << '\t' << cfg.contact_mode << '\t'
                   << cfg.protocol << '\t' << std::setprecision(17) << cfg.omega
                   << '\t' << cfg.delta << '\t' << cfg.beta_n << '\t' << cfg.beta_a << '\t' << cfg.beta_q
                   << '\t' << cfg.beta_mix_a << '\t' << cfg.beta_mixing_rule << '\t' << (cfg.claim_grade_corrected_model?1:0) << '\t' << cfg.tau << '\t' << cfg.i_c
                   << '\t' << cfg.fermi_k << '\t' << cfg.profile
                   << '\t' << cfg.seed_policy << '\t' << out.events << '\t'
                   << out.accepted_event_count << '\t' << int(out.censored)
                   << '\t' << out.stop_reason << '\t' << out.outbreak_size
                   << '\t' << out.ever_infected_count << '\t'
                   << out.duration << '\t' << out.peak_infected << '\t'
                   << out.peak_time << '\t' << out.max_q << '\t' << out.final_q
                   << '\t' << out.q_extinction_time << '\t' << out.q_person_time
                   << '\t' << out.q_interface_person_time
                   << '\t' << out.q_far_person_time
                   << '\t' << out.si_q_contact_time
                   << '\t' << out.si_n_contact_time
                   << '\t' << int(out.q_mesoscopic) << '\t' << int(out.q_macroscopic)
                   << '\t' << int(out.activated_ic) << '\t' << out.seed_vertex
                   << '\t' << out.seed_degree << '\t' << out.q_seed_vertex
                   << '\t' << out.q_seed_degree << '\t' << out.q_seed_distance
                   << '\t' << out.infected_degree_mean << '\t' << out.hubs_infected
                   << '\t' << cfg.n0 << '\t' << cfg.n_occ << '\t' << cfg.n_gc
                   << '\t' << cfg.behavior_mode << '\t' << cfg.q_placement
                   << '\t' << cfg.arm << '\t' << cfg.q0 << '\t' << cfg.i0
                   << '\t' << cfg.information_mode
                   << '\t' << cfg.information_radius
                   << '\t' << cfg.information_alpha
                   << '\t' << cfg.signal_seed
                   << '\t' << (cfg.instrumentation ? "on" : "off")
                   << '\t' << cfg.behavior_bin_width
                   << '\t' << out.initial_q << '\t' << out.initial_i
                   << '\t' << out.initial_iq << '\t' << out.min_q
                   << '\t' << out.first_ic_time << '\t' << out.q_at_first_ic
                   << '\t' << out.initial_q_lost_before_ic
                   << '\t' << out.time_above_ic << '\t' << out.ic_crossings
                   << '\t' << out.q_to_n_events << '\t' << out.n_to_q_events
                   << '\t' << out.infections_into_n << '\t' << out.infections_into_q
                   << '\t' << out.n_infection_hazard_time
                   << '\t' << out.q_infection_hazard_time
                   << '\t' << out.q_hazard_fraction
                   << '\t' << out.qn_local_prevalence_sum
                   << '\t' << out.nq_local_prevalence_sum
                   << '\t' << out.qn_near_infected_events
                   << '\t' << out.nq_near_infected_events
                   << '\t' << out.qn_mean_local_prevalence
                   << '\t' << out.nq_mean_local_prevalence
                   << '\t' << out.max_q_after_alert
                   << '\t' << out.q_recovery_after_alert
                   << '\t' << out.trajectory_signed_area
                   << '\t' << out.qn_max_bin << '\t' << out.nq_max_bin
                   << '\t' << out.qn_max_avalanche
                   << '\t' << out.nq_max_avalanche
                   << '\t' << out.qn_concentration
                   << '\t' << out.nq_concentration
                   << '\t' << out.qn_effective_bins
                   << '\t' << out.nq_effective_bins
                   << '\t' << out.qn_peak_time << '\t' << out.nq_peak_time
                   << '\t' << out.qn_peak_lag_from_i_peak
                   << '\t' << out.nq_peak_lag_from_i_peak
                   << '\t' << out.net_flux_peak
                   << '\t' << out.net_flux_peak_time
                   << '\t' << out.initial_q_components
                   << '\t' << out.initial_q_largest_component
                   << '\t' << out.initial_q_largest_fraction
                   << '\t' << out.initial_q_component_second_moment
                   << '\t' << out.initial_qn_interface
                   << '\t' << out.initial_qn_interface_fraction
                   << '\t' << out.seed_distance_to_q
                   << '\t' << out.seed_q_coverage_r1
                   << '\t' << out.seed_q_coverage_r2
                   << '\t' << out.seed_q_coverage_r4
                   << '\t' << out.null_behavior_proposals
                   << '\t' << out.max_ever_distance
                   << '\t' << out.signal_neighborhood_size
                   << '\t' << out.rms_ever_distance
                   << '\t' << (cfg.trace_mode ? "on" : "off")
                   << '\t' << out.qn_histogram << '\t' << out.nq_histogram
                   << '\t' << out.iq_trace
                   << '\t' << out.population_interface_trace
                   << '\t' << out.shell_trace << '\t' << out.trace_v2 << '\n';
        }
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
