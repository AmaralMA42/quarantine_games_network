#pragma once

#include <cstdint>
#include <limits>
#include <string>

namespace fpm {

// Claim-facing and diagnostic trajectory observables.
//
// Keep this structure semantically aligned with the output schema.  Adding an
// observer is allowed only if it is neutral to the stochastic dynamics; changing
// an existing definition requires a schema/version change and regression review.
struct Outcome {
    int64_t events = 0;
    int64_t accepted_event_count = 0;
    bool censored = false;
    int outbreak_size = 0;
    int ever_infected_count = 0;
    double duration = 0.0;
    int peak_infected = 0;
    double peak_time = 0.0;
    int max_q = 0;
    int final_q = 0;
    double q_extinction_time = std::numeric_limits<double>::quiet_NaN();
    double q_person_time = 0.0;
    double q_interface_person_time = 0.0;
    double q_far_person_time = 0.0;
    double si_q_contact_time = 0.0;
    double si_n_contact_time = 0.0;
    bool q_mesoscopic = false;
    bool q_macroscopic = false;
    bool activated_ic = false;
    int seed_vertex = -1;
    int seed_degree = 0;
    int q_seed_vertex = -1;
    int q_seed_degree = 0;
    int q_seed_distance = -1;
    double infected_degree_mean = std::numeric_limits<double>::quiet_NaN();
    int hubs_infected = 0;
    int initial_q = -1;
    int initial_i = -1;
    int initial_iq = -1;
    int min_q = -1;
    double first_ic_time = std::numeric_limits<double>::quiet_NaN();
    int q_at_first_ic = -1;
    double initial_q_lost_before_ic = std::numeric_limits<double>::quiet_NaN();
    double time_above_ic = 0.0;
    int ic_crossings = 0;
    int q_to_n_events = 0;
    int n_to_q_events = 0;
    int infections_into_n = 0;
    int infections_into_q = 0;
    double n_infection_hazard_time = 0.0;
    double q_infection_hazard_time = 0.0;
    double q_hazard_fraction = std::numeric_limits<double>::quiet_NaN();
    double qn_local_prevalence_sum = 0.0;
    double nq_local_prevalence_sum = 0.0;
    int qn_near_infected_events = 0;
    int nq_near_infected_events = 0;
    double qn_mean_local_prevalence = std::numeric_limits<double>::quiet_NaN();
    double nq_mean_local_prevalence = std::numeric_limits<double>::quiet_NaN();
    int max_q_after_alert = -1;
    double q_recovery_after_alert = std::numeric_limits<double>::quiet_NaN();
    double trajectory_signed_area = std::numeric_limits<double>::quiet_NaN();
    int qn_max_bin = 0;
    int nq_max_bin = 0;
    double qn_max_avalanche = std::numeric_limits<double>::quiet_NaN();
    double nq_max_avalanche = std::numeric_limits<double>::quiet_NaN();
    double qn_concentration = std::numeric_limits<double>::quiet_NaN();
    double nq_concentration = std::numeric_limits<double>::quiet_NaN();
    double qn_effective_bins = std::numeric_limits<double>::quiet_NaN();
    double nq_effective_bins = std::numeric_limits<double>::quiet_NaN();
    double qn_peak_time = std::numeric_limits<double>::quiet_NaN();
    double nq_peak_time = std::numeric_limits<double>::quiet_NaN();
    double qn_peak_lag_from_i_peak = std::numeric_limits<double>::quiet_NaN();
    double nq_peak_lag_from_i_peak = std::numeric_limits<double>::quiet_NaN();
    double net_flux_peak = std::numeric_limits<double>::quiet_NaN();
    double net_flux_peak_time = std::numeric_limits<double>::quiet_NaN();
    int initial_q_components = -1;
    int initial_q_largest_component = -1;
    double initial_q_largest_fraction = std::numeric_limits<double>::quiet_NaN();
    double initial_q_component_second_moment =
        std::numeric_limits<double>::quiet_NaN();
    int initial_qn_interface = -1;
    double initial_qn_interface_fraction = std::numeric_limits<double>::quiet_NaN();
    int seed_distance_to_q = -1;
    double seed_q_coverage_r1 = std::numeric_limits<double>::quiet_NaN();
    double seed_q_coverage_r2 = std::numeric_limits<double>::quiet_NaN();
    double seed_q_coverage_r4 = std::numeric_limits<double>::quiet_NaN();
    int null_behavior_proposals = 0;
    int signal_neighborhood_size = 0;
    int max_ever_distance = 0;
    double rms_ever_distance = 0.0;
    std::string qn_histogram;
    std::string nq_histogram;
    std::string iq_trace;
    std::string population_interface_trace;
    std::string shell_trace;
    std::string trace_v2;
    std::string stop_reason;
};

} // namespace fpm
