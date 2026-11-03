#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include "fpm/stable_rng.hpp"

namespace fpm {

/** @brief O(1) dense sparse-set pool for uniformly sampled integer items. */
class Pool {
public:
    std::vector<int> items;
    std::vector<int> position;
    explicit Pool(int capacity = 0) : position(capacity, -1) {}
    bool contains(int item) const { return position[item] >= 0; }
    void add(int item) {
        if (contains(item)) return;
        position[item] = int(items.size());
        items.push_back(item);
    }
    void remove(int item) {
        const int pos = position[item];
        if (pos < 0) return;
        const int last = items.back();
        items[pos] = last;
        position[last] = pos;
        items.pop_back();
        position[item] = -1;
    }
    int sample(StableRng &rng) const {
        return items[size_t(rng.bounded(items.size()))];
    }
    size_t size() const { return items.size(); }
};

/** @brief Reference weighted sampler preserving historical M20B event ordering. */
class Fenwick {
    std::vector<double> tree;
public:
    explicit Fenwick(int n = 0) : tree(n + 1, 0.0) {}
    void add(int index, double delta) {
        for (int i = index + 1; i < int(tree.size()); i += i & -i) tree[i] += delta;
    }
    double total() const {
        double result = 0.0;
        for (int i = int(tree.size()) - 1; i > 0; i -= i & -i) result += tree[i];
        return result;
    }
    int sample(double target) const {
        int index = 0;
        double accumulated = 0.0;
        int bit = 1;
        while ((bit << 1) < int(tree.size())) bit <<= 1;
        for (; bit; bit >>= 1) {
            const int next = index + bit;
            if (next < int(tree.size()) && accumulated + tree[next] <= target) {
                accumulated += tree[next];
                index = next;
            }
        }
        return std::min(index, int(tree.size()) - 2);
    }
};



/**
 * @brief Equal-weight sparse pool used only by the opt-in uniform_regular backend.
 *
 * It targets the same category-level stochastic law on compatible regular graphs,
 * but does not preserve historical seed-to-edge mapping.
 */
class UniformSparsePool {
    std::vector<int> items_;
    std::vector<int> position_;
    double weight_ = 0.0;
public:
    explicit UniformSparsePool(int capacity = 0, double weight = 0.0)
        : position_(capacity, -1), weight_(weight) {}
    void reset(int capacity, double weight) {
        items_.clear(); position_.assign(capacity, -1); weight_ = weight;
    }
    bool contains(int item) const { return position_[item] >= 0; }
    void add(int item) {
        if (contains(item)) return;
        position_[item] = int(items_.size()); items_.push_back(item);
    }
    void remove(int item) {
        const int pos = position_[item]; if (pos < 0) return;
        const int last = items_.back(); items_[pos] = last; position_[last] = pos;
        items_.pop_back(); position_[item] = -1;
    }
    double total() const { return weight_ * double(items_.size()); }
    int sample(double target) const {
        if (items_.empty() || !(weight_ > 0.0)) throw std::runtime_error("sample from empty uniform pool");
        size_t index = size_t(std::floor(target / weight_));
        if (index >= items_.size()) index = items_.size() - 1;
        return items_[index];
    }
};

/** @brief Edge sampling reproducibility regime. */
enum class EdgePoolBackend { FenwickExact, UniformRegular };

/** @brief Backend-neutral edge-category pool wrapper. */
class EdgePool {
    EdgePoolBackend backend_ = EdgePoolBackend::FenwickExact;
    Fenwick fenwick_;
    UniformSparsePool uniform_;
    double uniform_weight_ = 0.0;
public:
    EdgePool() = default;
    EdgePool(int capacity, EdgePoolBackend backend, double uniform_weight)
        : backend_(backend), uniform_weight_(uniform_weight) {
        if (backend_ == EdgePoolBackend::FenwickExact) fenwick_ = Fenwick(capacity);
        else uniform_.reset(capacity, uniform_weight);
    }
    void add(int item, double delta) {
        if (backend_ == EdgePoolBackend::FenwickExact) { fenwick_.add(item, delta); return; }
        if (delta == 0.0) return;
        const double magnitude = std::abs(delta);
        const double tolerance = 1e-13 * std::max(1.0, std::abs(uniform_weight_));
        if (std::abs(magnitude - uniform_weight_) > tolerance)
            throw std::runtime_error("uniform_regular received non-uniform edge weight");
        if (delta > 0.0) uniform_.add(item); else uniform_.remove(item);
    }
    double total() const { return backend_ == EdgePoolBackend::FenwickExact ? fenwick_.total() : uniform_.total(); }
    int sample(double target) const { return backend_ == EdgePoolBackend::FenwickExact ? fenwick_.sample(target) : uniform_.sample(target); }
};

} // namespace fpm
