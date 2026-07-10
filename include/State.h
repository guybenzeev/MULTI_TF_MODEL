#ifndef STATE_H
#define STATE_H

#include "Edit.h"
#include <functional>
#include <cstddef>

enum class StateType {
    FREE,
    BOUND_NS,
    BOUND_S
};

/**
 * @brief Represents the state of one sub-chain site.
 *
 * protein stores the index + 1 of the protein in SubChainTopo::proteins_.
 * protein == 0 means no protein is bound.
 */
struct State {
    StateType state;
    int strandSide;
    int protein;

    State(StateType state_val = StateType::FREE, int strand_side_val = 0, int protein_val = 0)
        : state(state_val), strandSide(strand_side_val), protein(protein_val)
    {}

    void changeState(const Edit& edit) {
        switch(edit.type) {
            case EditType::BIND_NS:
                bind_ns(edit.strandSide, edit.protein);
                break;
            case EditType::BIND_S:
                bind_s();
                break;
            case EditType::UNBIND_NS:
                unbind();
                break;
            case EditType::UNBIND_S:
                unbind_s();
                break;
            case EditType::SWITCH_SIDE:
                switch_side(edit.strandSide);
                break;
            case EditType::SLIDE_RIGHT:
            case EditType::SLIDE_LEFT:
                slideFrom();
                break;
            default:
                break;
        }
    }

    bool isFree() const {
        return state == StateType::FREE;
    }

    bool isBoundNs() const {
        return state == StateType::BOUND_NS;
    }

    bool isBoundS() const {
        return state == StateType::BOUND_S;
    }

    void unbind() {
        if (!isBoundNs()) {
            return;
        }
        state = StateType::FREE;
        strandSide = 0;
        protein = 0;
    }

    void bind_ns(int side, int protein_index) {
        if (!isFree() || side == 0 || protein_index == 0) {
            return;
        }
        state = StateType::BOUND_NS;
        strandSide = side;
        protein = protein_index;
    }

    void bind_s() {
        if (!isBoundNs()) {
            return;
        }
        state = StateType::BOUND_S;
    }

    void unbind_s() {
        if (!isBoundS()) {
            return;
        }
        state = StateType::BOUND_NS;
    }

    void switch_side(int side) {
        if (!isBoundNs() || side == 0) {
            return;
        }
        strandSide = side;
    }

    void slideFrom() {
        if (!isBoundNs()) {
            return;
        }
        state = StateType::FREE;
        strandSide = 0;
        protein = 0;
    }

    void slideTo(int side, int protein_index) {
        bind_ns(side, protein_index);
    }
    
    int getStateID(int num_sides, int num_proteins) const {
        if (isFree()) {
            return 0;
        }

        if (protein <= 0 || protein > num_proteins ||
            strandSide <= 0 || strandSide > num_sides) {
            return -1;
        }

        int proteinIndex = protein - 1;
        int sideIndex = strandSide - 1;

        int perBindingType = num_sides * num_proteins;

        if (isBoundNs()) {
            return 1 + proteinIndex * num_sides + sideIndex;
        }

        if (isBoundS()) {
            return 1 + perBindingType + proteinIndex * num_sides + sideIndex;
        }

        return -1;
    }

    bool operator==(const State& other) const {
        return state == other.state &&
            strandSide == other.strandSide &&
            protein == other.protein;
    }

};


struct StateHash {
    std::size_t operator()(const State& s) const {
        std::size_t h1 = std::hash<int>{}(static_cast<int>(s.state));
        std::size_t h2 = std::hash<int>{}(s.strandSide);
        std::size_t h3 = std::hash<int>{}(s.protein);

        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

#endif
