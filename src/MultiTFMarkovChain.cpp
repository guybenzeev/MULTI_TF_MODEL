// MultiTFMarkovChain.cpp
#include "MultiTFMarkovChain.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <iomanip>


// ---------- Private helpers (stubs) ----------


double MultiTFMarkovChain::findHoldingTime() {
    if (totalExitRate_ <= 0.0) {
        throw std::runtime_error("Total exit rate is zero or negative, cannot sample holding time.");
    }

    std::exponential_distribution<double> dist(totalExitRate_);
    return dist(rng_);
}

int MultiTFMarkovChain::selectSubChain(double& threshold) {
    double cumulativeRate = 0.0;
    double prev;

    for(int i = 0; i < length_; ++i) {
        prev = cumulativeRate;
        cumulativeRate += subChainExitRates_[i];
        if (threshold < cumulativeRate) {
            threshold -= prev; // update threshold to the start of this interval
            return i;
        }
    }
    return -1;
}

int MultiTFMarkovChain::updatePossibleEdits() {
    // Stub: no rate updates yet
    totalExitRate_ = 0.0;
    subChainExitRates_.clear();
    return 0;
}

void MultiTFMarkovChain::updateChain() {
    // Stub: no simulation step yet
}

void MultiTFMarkovChain::initializeExitRates(){
    subChainExitRates_.clear();
    subChainExitRates_.reserve(length_);
    totalExitRate_ = 0.0;
    for(int i = 0; i < length_; ++i){
        double rate = currentChain_[i]->getTotalExitRate(currentChain_);
        subChainExitRates_.push_back(rate);
        totalExitRate_ += rate;
    }
}

void MultiTFMarkovChain::updateExitRates(int centerIndex){
    int start = std::max(0, centerIndex - effectRadius_);
    int end = std::min(length_ - 1, centerIndex + effectRadius_);

    for(int i = start; i <= end; ++i){
        double oldRate = subChainExitRates_[i];
        double newRate = currentChain_[i]->getTotalExitRate(currentChain_);
        subChainExitRates_[i] = newRate;
        totalExitRate_ += (newRate - oldRate);
    }
}




// ---------- Constructors / Destructor ----------

MultiTFMarkovChain::MultiTFMarkovChain(int length, const SubChainTopo& topo)
    : topo_(topo),
      length_(length),
      currentTime_(0),
      pastEvents_(),
      pastTimes_(),
      effectRadius_(topo.getEffectRadius()),
      subChainExitRates_(),
      totalExitRate_(0.0),
      rng_(std::random_device{}()) 
      {
        currentChain_.reserve(length_);

        for (int i = 0; i < length_; ++i) {
            currentChain_.push_back(
                std::make_unique<SubChain>(
                    topo_,
                    i,
                    0  // initial state ID = 0 (free)
                )
            );
        }

        initializeExitRates();
      }

MultiTFMarkovChain::~MultiTFMarkovChain() = default;

// ---------- Public API ----------

void MultiTFMarkovChain::runSimulation(double runTime) {
    pastEvents_.clear();
    pastTimes_.clear();
    
    currentTime_ = 0;
    initializeExitRates();
    Edit slidToEdit;
    slidToEdit.type = EditType::BIND_NS;

    //Clear nextEdit?
    while (currentTime_ < runTime) {
        std::cout << "Total exit rate: " << totalExitRate_ << "\n";
        printCurrentStates();
        double holdingTime = findHoldingTime(); // sample holding time
        currentTime_ += holdingTime;

        if (currentTime_ > runTime) break;

        // Select which sub-chain will fire next
        std::uniform_real_distribution<double> dist(0.0, totalExitRate_);
        double threshold = dist(rng_);
        selectedSubChainIndex_ = selectSubChain(threshold);

        if (selectedSubChainIndex_ < 0) {
            throw std::runtime_error("selectSubChain failed (threshold out of range?)");
        }

        // Find and apply the next edit
        Edit nextEdit = currentChain_[selectedSubChainIndex_]->findNextEdit(threshold, currentChain_);
        if(nextEdit.type == EditType::SLIDE_LEFT || nextEdit.type == EditType::SLIDE_RIGHT){
            std::cout << "Selected edit: " << nextEdit.toString() << " on sub-chain " << selectedSubChainIndex_ << "\n";
            printCurrentStates();
        }
        currentChain_[selectedSubChainIndex_]->applyEdit(nextEdit);

        //CHECK THIS PART FOR SLIDING UPDATES
        if(nextEdit.type == EditType::SLIDE_LEFT) {
            slidToEdit.strandSide = nextEdit.strandSide;
            currentChain_[selectedSubChainIndex_ - 1]->applyEdit(slidToEdit);
            updateExitRates(selectedSubChainIndex_ - 1);
        }
        else if(nextEdit.type == EditType::SLIDE_RIGHT) {
            slidToEdit.strandSide = nextEdit.strandSide;
            currentChain_[selectedSubChainIndex_ + 1]->applyEdit(slidToEdit);
            updateExitRates(selectedSubChainIndex_ + 1);
        }
        if(nextEdit.type == EditType::SLIDE_LEFT || nextEdit.type == EditType::SLIDE_RIGHT){
            std::cout << "after edit: " << selectedSubChainIndex_ << "\n";
            printCurrentStates();
        }

        // Update exit rates
        updateExitRates(selectedSubChainIndex_);

        // Record event time + event
        pastTimes_.push_back(currentTime_);
        pastEvents_.push_back(PastEvent{selectedSubChainIndex_, nextEdit});
    }
}

void MultiTFMarkovChain::printCurrentStates() const{
    for (const auto& subChainPtr : currentChain_) {
        const State& state = subChainPtr->getCurrentState();
        int stateID = state.getStateID(topo_.getNumSides());
        std::cout << stateID << ",";
    }
    std::cout << "\n";
}

void MultiTFMarkovChain::printExitRates(){
    initializeExitRates();
    std::cout << "Sub-chain Exit Rates:\n";
    for (size_t i = 0; i < subChainExitRates_.size(); ++i) {
        std::cout << "  Sub-chain " << i << ": " << subChainExitRates_[i] << "\n";
    }
    std::cout << "Total Exit Rate: " << totalExitRate_ << "\n";
}

double MultiTFMarkovChain::getTotalExitRate() const {
    return totalExitRate_;
}

const std::vector<double>& MultiTFMarkovChain::getSubChainExitRates() const {
    return subChainExitRates_;
}

void MultiTFMarkovChain::writeToCSV(const std::string& filename, double runTime) const {
    std::ofstream out(filename);

    if (!out) {
        throw std::runtime_error("Failed to open CSV file: " + filename);
    }

    // Header
   // out << "time,index,edit_type,side\n";
   out << "time,index,edit,side,length=" << length_
    << ",num_sides=" << topo_.getNumSides() << ",runTime=" << runTime << "\n";

    for (size_t i = 0; i < pastTimes_.size(); ++i) {
        const auto& ev = pastEvents_[i];

        out << std::fixed << std::setprecision(10) << pastTimes_[i] << ","
            << ev.index << ","
            //<< static_cast<int>(ev.edit.type) << ","
            << ev.edit.typeAsString() << ","
            << ev.edit.strandSide << "\n";
    }

    out.close();
}

