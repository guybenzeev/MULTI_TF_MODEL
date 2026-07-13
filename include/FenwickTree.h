#pragma once

#include <vector>

class FenwickTree {
private:
    std::vector<double> tree_;

public:
    explicit FenwickTree(int size);

    void add(int index, double difference);
    double prefixSum(int index) const;
    double totalSum() const;
    int findIndex(double threshold) const;
};