#ifndef STATE_H
#define STATE_H

#include <cstddef>
#include <functional>

#include "Edit.h"

/**
 * @file State.h
 * @brief Declares site states and their hash function.
 */

/**
 * @enum StateType
 * @brief Identifies the binding state at a strand site.
 */
enum class StateType {
    FREE,     /**< No protein is bound at the site. */
    BOUND_NS, /**< A protein is bound non-specifically. */
    BOUND_S   /**< A protein is bound specifically. */
};

/**
 * @struct State
 * @brief Represents the state of one sub-chain site.
 *
 * `protein` stores a one-based protein identifier. A value of zero means that
 * no protein is bound. Multi-site proteins mark their non-head sites occupied.
 */
struct State {
    StateType state; /**< Current binding category. */
    int strandSide;  /**< One-based strand side, or zero when free. */
    int protein;     /**< One-based protein identifier, or zero when free. */
    bool occupied;   /**< Whether this is a non-head site in a protein footprint. */

    /**
     * @brief Constructs a site state.
     * @param state_val Initial binding category.
     * @param strand_side_val Initial one-based strand side.
     * @param protein_val Initial one-based protein identifier.
     */
    State(StateType state_val = StateType::FREE, int strand_side_val = 0, int protein_val = 0)
        : state(state_val), strandSide(strand_side_val), protein(protein_val), occupied(false)
    {}

    /**
     * @brief Applies an edit's local state transition.
     * @param edit Transition to apply.
     */
    void changeState(const Edit& edit) {
        occupied = false;
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

    /** @return `true` when the site is free. */
    bool isFree() const {
        return state == StateType::FREE;
    }

    /** @return `true` when the site is bound non-specifically. */
    bool isBoundNs() const {
        return state == StateType::BOUND_NS;
    }

    /** @return `true` when the site is bound specifically. */
    bool isBoundS() const {
        return state == StateType::BOUND_S;
    }

    /** @return `true` for a non-head site occupied by a multi-site protein. */
    bool isOccupied() const {
        return occupied;
    }

    /** @brief Marks the site as a non-head footprint site. */
    void occupy() {
        occupied = true;
    }

    /** @brief Releases a non-specifically bound protein. */
    void unbind() {
        if (!isBoundNs()) {
            return;
        }
        state = StateType::FREE;
        strandSide = 0;
        protein = 0;
    }

    /**
     * @brief Binds a protein non-specifically when the site is free.
     * @param side One-based strand side.
     * @param protein_index One-based protein identifier.
     */
    void bind_ns(int side, int protein_index) {
        if (!isFree() || side == 0 || protein_index == 0) {
            return;
        }
        state = StateType::BOUND_NS;
        strandSide = side;
        protein = protein_index;
    }

    /** @brief Converts a non-specific binding to a specific binding. */
    void bind_s() {
        if (!isBoundNs()) {
            return;
        }
        state = StateType::BOUND_S;
    }

    /** @brief Converts a specific binding to a non-specific binding. */
    void unbind_s() {
        if (!isBoundS()) {
            return;
        }
        state = StateType::BOUND_NS;
    }

    /**
     * @brief Changes the side of a non-specifically bound protein.
     * @param side New one-based strand side.
     */
    void switch_side(int side) {
        if (!isBoundNs() || side == 0) {
            return;
        }
        strandSide = side;
    }

    /** @brief Clears the source site of a sliding protein. */
    void slideFrom() {
        if (!isBoundNs()) {
            return;
        }
        state = StateType::FREE;
        strandSide = 0;
        protein = 0;
    }

    /**
     * @brief Converts this state into its compact topology-specific identifier.
     * @param num_sides Number of strand sides in the topology.
     * @param num_proteins Number of protein types in the topology.
     * @return Zero for a free site, a positive state ID, or -1 for invalid data.
     */
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

    /**
     * @brief Compares the binding identity of two states.
     *
     * The occupied marker is intentionally excluded because edit tables are
     * keyed by binding state, side, and protein.
     *
     * @param other State to compare.
     * @return `true` when the binding fields are equal.
     */
    bool operator==(const State& other) const {
        return state == other.state &&
            strandSide == other.strandSide &&
            protein == other.protein;
    }
};

/**
 * @struct StateHash
 * @brief Hashes a State for use in topology edit dictionaries.
 */
struct StateHash {
    /**
     * @brief Computes a hash from a state's binding identity.
     * @param s State to hash.
     * @return Hash value compatible with State::operator==.
     */
    std::size_t operator()(const State& s) const {
        std::size_t h1 = std::hash<int>{}(static_cast<int>(s.state));
        std::size_t h2 = std::hash<int>{}(s.strandSide);
        std::size_t h3 = std::hash<int>{}(s.protein);

        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

#endif // STATE_H
