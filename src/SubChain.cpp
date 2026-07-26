#include "SubChain.h"
#include <stdexcept>
#include <iostream>

Edit SubChain::findNextEdit(double threshold, const std::vector<std::unique_ptr<SubChain>>& chain) {
    const auto& edits = getPossibleEditsForCurrentState();
    if (edits.empty()) {
        throw std::runtime_error("No possible edits for current state");
    }

    const State& current = currentState;

    double cumulativeRate = 0.0;
    double rate = 0.0;
    
    //std::cout << "selected index: " << nodeId_
    //          << ", current State: "
    //          << currentState.getStateID(topo_.getNumSides(), topo_.getNumProteins())
    //          << "\n";


    for (const Edit& e : edits) {
        rate = topo_.computeRate(current, e, chain, nodeId_);
        if (rate > 0.0) {
            cumulativeRate += rate;
        }
        
        //std::cout << "Edit type: " << e.toString() << ", rate: " << rate << ", cumulativeRate: " << cumulativeRate <<", Threshold: " << threshold << "\n";

        if (cumulativeRate >= threshold) {
            return e;
        }
    }
    throw std::runtime_error("reached end of findNextEdit without selecting an edit, incorrect threshold?");
}

double SubChain::getTotalExitRate(const std::vector<std::unique_ptr<SubChain>>& chain) {
    if(currentState.isOccupied()) {
        return 0.0;
    }
    
    const auto& edits = getPossibleEditsForCurrentState();
    if (edits.empty()) {
        return 0.0;
    }

    const State& current = currentState;

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

void SubChain::applyEdit(Edit& edit, int proteinHeadIndex) {
    currentState.changeState(edit);
    if (nodeId_ != proteinHeadIndex) {
        currentState.occupy();
    }
}
