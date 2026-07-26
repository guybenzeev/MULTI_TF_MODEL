#pragma once

#include <vector>

/**
 * @file FenwickTree.h
 * @brief Declares a Fenwick tree for cumulative transition rates.
 */

/**
 * @class FenwickTree
 * @brief Maintains prefix sums while supporting logarithmic point updates.
 *
 * Public indices are zero-based. The tree stores non-negative transition
 * rates and is used to select a sub-chain from a cumulative-rate threshold.
 */
class FenwickTree {
private:
    std::vector<double> tree_; /**< One-based internal Fenwick-tree storage. */

public:
    /**
     * @brief Creates an empty rate tree.
     * @param size Number of values represented by the tree.
     * @pre @p size is non-negative.
     */
    explicit FenwickTree(int size);

    /**
     * @brief Adds a difference to one value.
     * @param index Zero-based value index.
     * @param difference Signed amount to add.
     */
    void add(int index, double difference);

    /**
     * @brief Computes the inclusive prefix sum.
     * @param index Zero-based last index in the prefix.
     * @return Sum of values in the range `[0, index]`.
     */
    double prefixSum(int index) const;

    /**
     * @brief Locates the cumulative-rate interval containing a threshold.
     * @param threshold Zero-based cumulative threshold.
     * @return Zero-based index whose interval contains @p threshold.
     */
    int findIndex(double threshold) const;
};
