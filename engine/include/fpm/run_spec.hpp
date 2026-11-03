#pragma once

#include <cstdint>
#include <limits>
#include <string>

namespace fpm {

// Fully materialized per-trajectory run specification.
//
// The parser below is deliberately strict: a trajectory is scientifically
// meaningful only together with its schema, graph identity and domain-separated
// random streams.  Mission names/stages are metadata; they never select a
// different simulator implementation.
struct Config {
    int run_index{};
    std::string schema_version, analysis_role, run_fingerprint, fingerprint_version;
    std::string run_id, stage, topology, graph_path, graph_sha256;
    std::string graph_certificate_path, graph_certificate_sha256, boundary_condition;
    std::string contact_mode, profile, arm;
    std::string protocol = "LEGACY_IC_FIXED";
    int population{};
    int n0{}, n_occ{}, n_gc{};
    uint64_t graph_seed{}, q_init_seed{}, i_init_seed{}, dynamic_seed{};
    double beta_n{}, beta_a{}, beta_q{}, beta_mix_a{}, q0{}, i0{}, tau{}, i_c{}, fermi_k{};
    std::string beta_mixing_rule = "INTERPOLATED";
    bool claim_grade_corrected_model = true;
    bool claim_grade_materiality = true;
    double omega = 1.0;
    double delta = std::numeric_limits<double>::quiet_NaN();
    std::string seed_policy, behavior_mode, q_placement;
    std::string information_mode = "global";
    int information_radius = 1;
    double information_alpha = 0.0;
    uint64_t signal_seed{};
    bool instrumentation = true;
    bool trace_mode = false;
    double behavior_bin_width = 0.25;
    int64_t event_limit{};
    double tmax{};
};

} // namespace fpm
