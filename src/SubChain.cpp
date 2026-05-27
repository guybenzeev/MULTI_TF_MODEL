#include "SubChain.h"
#include <stdexcept>
#include <iostream>


double SubChain::computeInternalRate(
    int,
    const std::vector<std::unique_ptr<SubChain>>&
) {
    return 0.0;
}

double SubChain::computeSlidingRate(
    int,
    const std::vector<std::unique_ptr<SubChain>>&
) {
    return 0.0;
}

Edit SubChain::findNextEdit(double threshold, const std::vector<std::unique_ptr<State>>& chain) {
    const auto& edits = getPossibleEditsForCurrentState();
    if (edits.empty()) {
        throw std::runtime_error("No possible edits for current state");
    }

    const State& current = *chain.at(static_cast<size_t>(nodeId_));

    double cumulativeRate = 0.0;
    double rate = 0.0;


    for (const Edit& e : edits) {
        rate = topo_.computeRate(current, e, chain, nodeId_);
        if (rate > 0.0) {
            cumulativeRate += rate;
        }
        std::cout << "Edit type: " << e.toString() << ", rate: " << rate << ", cumulativeRate: " << cumulativeRate <<", Threshold: " << threshold << "\n";

        if (cumulativeRate >= threshold) {
            return e;
        }
    }
    throw std::runtime_error("reached end of findNextEdit without selecting an edit, incorrect threshold?");
}

double SubChain::getTotalExitRate(const std::vector<std::unique_ptr<State>>& chain) {
    // currentStateID is your "stateId" index into the edit table
    auto it = editTable_.find(currentStateID);
    if (it == editTable_.end()) {
        return 0.0;
    }

    const auto& edits = it->second;
    if (edits.empty()) {
        return 0.0;
    }

    // Current local state for this node is the State stored in the global chain at nodeId_
    // (This matches the topo::computeRate signature and lets rates depend on neighbors.)
    const State& current = *chain.at(static_cast<size_t>(nodeId_));

    double total = 0.0;
    for (const Edit& e : edits) {
        // topo_ decides which rate function to call based on e.type
        double r = topo_.computeRate(current, e, chain, nodeId_);
        if (r > 0.0) {
            total += r;
        }
    }

    return total;
}

void SubChain::applyEdit(Edit& edit) {
    currentState.changeState(edit); // Update internal State
    currentStateID = currentState.getStateID(topo_.getNumSides());
}
