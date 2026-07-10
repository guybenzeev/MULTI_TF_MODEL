#include "SingleStrandTopo.h"
#include "SubChain.h"

void SingleStrandTopo::buildPossibleEditsDict() {
    possibleEditsByState_.clear();

    // FREE state: can bind non-specific on any side with any protein.
    {
        State freeState{};
        possibleEditsByState_[freeState] = EditList{};

        for (int protein = 1; protein <= getNumProteins(); ++protein) {
            Edit e{};
            e.type = EditType::BIND_NS;
            e.strandSide = 1;
            e.protein = protein;
            possibleEditsByState_[freeState].push_back(e);
        }
    }

    // NS states: can slide, bind specific, unbind NS, switch side
    for (int protein = 1; protein <= getNumProteins(); ++protein) {
        State nsState{StateType::BOUND_NS, 1, protein};
        possibleEditsByState_[nsState] = EditList{};

        // Slide left
        {
            Edit e{};
            e.type = EditType::SLIDE_LEFT;
            e.strandSide = 1;
            e.protein = protein;
            possibleEditsByState_[nsState].push_back(e);
        }

        // Slide right
        {
            Edit e{};
            e.type = EditType::SLIDE_RIGHT;
            e.strandSide = 1;
            e.protein = protein;
            possibleEditsByState_[nsState].push_back(e);
        }

        // Bind specific
        {
            Edit e{};
            e.type = EditType::BIND_S;
            e.strandSide = 1;
            e.protein = protein;
            possibleEditsByState_[nsState].push_back(e);
        }

        // Unbind non-specific
        {
            Edit e{};
            e.type = EditType::UNBIND_NS;
            e.strandSide = 1;
            e.protein = protein;
            possibleEditsByState_[nsState].push_back(e);
        }
            
    }

    // S states: can unbind S (back to NS)
    for (int protein = 1; protein <= getNumProteins(); ++protein) {
        State sState{StateType::BOUND_S, 1, protein};
        possibleEditsByState_[sState] = EditList{};

        Edit e{};
        e.type = EditType::UNBIND_S;
        e.strandSide = 1;
        e.protein = protein;
        possibleEditsByState_[sState].push_back(e);
    }
}



double SingleStrandTopo::switchSideRate(
    const State& currentState,
    const Edit& edit,
    const std::vector<std::unique_ptr<SubChain>>& chain,
    const int nodeId
) {
    return 0.0;
}

double SingleStrandTopo::computeRate(
        const State& currentState,
        const Edit& edit,
        const std::vector<std::unique_ptr<SubChain>>& chain,
        const int nodeId

    ) const {
        if(edit.type == EditType::SWITCH_SIDE){
            return 0.0;
        }

        double rate = proteins_[edit.protein - 1]->computeRate(currentState, edit, chain, nodeId);
        return rate;
    }

