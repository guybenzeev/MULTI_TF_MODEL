#ifndef STATE_H
#define STATE_H

#include <string>
#include "Edit.h"

/**
 * @brief Represents a state in the Markov chain.
 *
 * A State is identified by a human-readable name and a numeric ID.
 */

struct State {
    int free;
    int ns;
    int s;

    State(int free_val = 1, int ns_val = 0, int s_val = 0)
        : free(free_val), ns(ns_val), s(s_val)
    {}

    void changeState(const Edit& edit){
        switch(edit.type){
            case EditType::BIND_NS:
                bind_ns(edit.strandSide);
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


    void unbind(){
        if(free != 0 || ns == 0 || s != 0){
            //TODO: implement error handling for incorrect state transition
            return;
        }
        ns = 0;
        free = 1;
    }

    void bind_ns(int side){
        if(free == 0 || ns != 0 || s != 0){
            //TODO: implement error handling for incorrect state transition
            return;
        }
        free = 0;
        ns = side;
    }

    void bind_s(){
        if(free != 0 || ns == 0 || s != 0){
            //TODO: implement error handling for incorrect state transition
            return;
        }
        s = ns;
        ns = 0;
    }

    void unbind_s(){
        if(s == 0 || free != 0 || ns != 0){
            //TODO: implement error handling for incorrect state transition
            return;
        }
        ns = s;
        s = 0;
    }

    void switch_side(int side){
        if(ns == 0 || free != 0 || s != 0){
            //TODO: implement error handling for incorrect state transition
            return;
        }
        ns = side;
    }

    void slideFrom(){
        if(s != 0 || free != 0 || ns == 0){
            //TODO: implement error handling for incorrect state transition
            return;
        }
        ns = 0;
        free = 1;
    }

    void slideTo(int side){
        bind_ns(side);
    }
    
    int getStateID(int num_sides) const{
        if(free == 1 && ns == 0 && s == 0){
            return 0;
        }
        else if(free == 0 && ns != 0 && s == 0){
            return ns;
        }
        else if(free == 0 && ns == 0 && s != 0){
            return num_sides + s;
        }
        else{
            return -1;
        }
    }
};



#endif