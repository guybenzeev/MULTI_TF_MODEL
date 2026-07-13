#include "FenwickTree.h"

FenwickTree::FenwickTree(int size)
    : tree_(size + 1, 0.0) {
}

void FenwickTree::add(int index, double change) {
    ++index;

    while (index < static_cast<int>(tree_.size())) {
        tree_[index] += change;
        index += index & -index;
    }
}

double FenwickTree::prefixSum(int index) const {
    double sum = 0.0;
    ++index;

    while (index > 0) {
        sum += tree_[index];
        index -= index & -index;
    }

    return sum;
}

double FenwickTree::totalSum() const {
    return prefixSum(
        static_cast<int>(tree_.size()) - 2
    );
}

int FenwickTree::findIndex(double threshold) const {
    int index = 0;
    double accumulated = 0.0;

    int bit = 1;

    while (bit < static_cast<int>(tree_.size())) {
        bit <<= 1;
    }

    for (bit >>= 1; bit > 0; bit >>= 1) {
        int next = index + bit;

        if (
            next < static_cast<int>(tree_.size()) &&
            accumulated + tree_[next] <= threshold
        ) {
            index = next;
            accumulated += tree_[next];
        }
    }

    return index;
}
