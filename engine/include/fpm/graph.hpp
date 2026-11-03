#pragma once

#include <algorithm>
#include <utility>
#include <vector>

namespace fpm {

/** @brief Validated simple undirected graph and lightweight structural metadata. */
struct Graph {
    int n = 0;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::vector<int>> neighbors;
    std::vector<int> degree;
    std::vector<char> hub;
};

/**
 * @brief Immutable graph-derived indexing shared by all trajectories on one graph.
 *
 * Directed-edge order intentionally matches the corrected M20B reference kernel,
 * allowing exact-backend refactors to preserve the historical seed mapping.
 */
struct GraphRuntime {
    std::vector<int> directed_u;
    std::vector<int> directed_v;
    std::vector<int> incident_offsets;
    std::vector<int> incident_edges;
    bool regular = false;
    int regular_degree = -1;
};

/** Build directed-edge and CSR-like incident indexing once per loaded graph. */
inline GraphRuntime build_graph_runtime(const Graph &graph) {
    GraphRuntime runtime;
    const int directed = 2 * int(graph.edges.size());
    runtime.directed_u.resize(directed);
    runtime.directed_v.resize(directed);
    runtime.incident_offsets.assign(graph.n + 1, 0);
    for (int u = 0; u < graph.n; ++u)
        runtime.incident_offsets[u + 1] = runtime.incident_offsets[u] + graph.degree[u];
    runtime.incident_edges.resize(directed);
    std::vector<int> cursor = runtime.incident_offsets;
    for (int e = 0; e < int(graph.edges.size()); ++e) {
        const auto [u, v] = graph.edges[e];
        runtime.directed_u[2*e] = u; runtime.directed_v[2*e] = v;
        runtime.directed_u[2*e+1] = v; runtime.directed_v[2*e+1] = u;
        runtime.incident_edges[cursor[u]++] = 2*e;
        runtime.incident_edges[cursor[v]++] = 2*e+1;
    }
    runtime.regular = !graph.degree.empty() && std::all_of(
        graph.degree.begin(), graph.degree.end(),
        [&](int d){ return d == graph.degree.front(); });
    runtime.regular_degree = runtime.regular ? graph.degree.front() : -1;
    return runtime;
}

} // namespace fpm
