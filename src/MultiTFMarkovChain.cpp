// MultiTFMarkovChain.cpp
#include "MultiTFMarkovChain.h"
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <iomanip>


namespace {
class CoutRedirect {
public:
    explicit CoutRedirect(std::ostream& target)
        : oldBuffer_(std::cout.rdbuf(target.rdbuf()))
    {}

    ~CoutRedirect() {
        std::cout.rdbuf(oldBuffer_);
    }

private:
    std::streambuf* oldBuffer_;
};
}


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
    for(int i = 0; i < length_; ++i){
        int distance = (i > centerIndex) ? (i - centerIndex) : (centerIndex - i);
        if (distance > currentChain_[i]->getRadius()) {
            continue;
        }

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
    std::ofstream coutLog("data/multitf_cout.log");
    if (!coutLog) {
        throw std::runtime_error("Failed to open cout log file: data/multitf_cout.log");
    }
    CoutRedirect coutRedirect(coutLog);

    pastEvents_.clear();
    pastTimes_.clear();
    
    currentTime_ = 0;
    initializeExitRates();
    Edit slidToEdit;
    slidToEdit.type = EditType::BIND_NS;
    int step = 0;

    //Clear nextEdit?
    while (currentTime_ < runTime) {
        std::cout << "\nSTEP " << step << "\n";
        std::cout << "TIME before=" << currentTime_ << "\n";
        std::cout << "RATES total_exit_rate=" << totalExitRate_ << "\n";
        std::cout << "RATES subchain_exit_rates=";
        for (double rate : subChainExitRates_) {
            std::cout << rate << ",";
        }
        std::cout << "\n";
        std::cout << "STATE before=";
        printCurrentStates();

        double holdingTime = findHoldingTime(); // sample holding time
        std::cout << "TIME holding_time=" << holdingTime << "\n";
        currentTime_ += holdingTime;
        std::cout << "TIME after_holding=" << currentTime_ << "\n";

        if (currentTime_ > runTime) {
            std::cout << "STEP_RESULT exceeded_run_time runTime=" << runTime << "\n";
            break;
        }

        // Select which sub-chain will fire next
        std::uniform_real_distribution<double> dist(0.0, totalExitRate_);
        double threshold = dist(rng_);
        double rawThreshold = threshold;
        std::cout << "THRESHOLD raw=" << rawThreshold << "\n";
        selectedSubChainIndex_ = selectSubChain(threshold);
        std::cout << "THRESHOLD local_after_subchain_selection=" << threshold << "\n";

        if (selectedSubChainIndex_ < 0) {
            throw std::runtime_error("selectSubChain failed (threshold out of range?)");
        }
        std::cout << "SELECTED subchain_index=" << selectedSubChainIndex_
                  << ", subchain_exit_rate=" << subChainExitRates_[selectedSubChainIndex_]
                  << "\n";

        // Find and apply the next edit
        Edit nextEdit = currentChain_[selectedSubChainIndex_]->findNextEdit(threshold, currentChain_);
        std::cout << "EDIT selected=" << nextEdit.toString()
                  << ", protein=" << nextEdit.protein
                  << "\n";
        if(nextEdit.type == EditType::SLIDE_LEFT || nextEdit.type == EditType::SLIDE_RIGHT){
            std::cout << "STATE before_slide=";
            printCurrentStates();
        }
        int slidingProtein = currentChain_[selectedSubChainIndex_]->getCurrentState().protein;
        std::cout << "STATE selected_subchain_before_apply="
                  << currentChain_[selectedSubChainIndex_]->getCurrentState().getStateID(topo_.getNumSides(), topo_.getNumProteins())
                  << ", sliding_protein=" << slidingProtein
                  << "\n";
        currentChain_[selectedSubChainIndex_]->applyEdit(nextEdit);
        std::cout << "STATE selected_subchain_after_apply="
                  << currentChain_[selectedSubChainIndex_]->getCurrentState().getStateID(topo_.getNumSides(), topo_.getNumProteins())
                  << "\n";

        //CHECK THIS PART FOR SLIDING UPDATES
        if(nextEdit.type == EditType::SLIDE_LEFT) {
            slidToEdit.strandSide = nextEdit.strandSide;
            slidToEdit.protein = slidingProtein;
            currentChain_[selectedSubChainIndex_ - 1]->applyEdit(slidToEdit);
            std::cout << "SLIDE destination_index=" << (selectedSubChainIndex_ - 1)
                      << ", destination_state="
                      << currentChain_[selectedSubChainIndex_ - 1]->getCurrentState().getStateID(topo_.getNumSides(), topo_.getNumProteins())
                      << "\n";
            updateExitRates(selectedSubChainIndex_ - 1);
        }
        else if(nextEdit.type == EditType::SLIDE_RIGHT) {
            slidToEdit.strandSide = nextEdit.strandSide;
            slidToEdit.protein = slidingProtein;
            currentChain_[selectedSubChainIndex_ + 1]->applyEdit(slidToEdit);
            std::cout << "SLIDE destination_index=" << (selectedSubChainIndex_ + 1)
                      << ", destination_state="
                      << currentChain_[selectedSubChainIndex_ + 1]->getCurrentState().getStateID(topo_.getNumSides(), topo_.getNumProteins())
                      << "\n";
            updateExitRates(selectedSubChainIndex_ + 1);
        }
        if(nextEdit.type == EditType::SLIDE_LEFT || nextEdit.type == EditType::SLIDE_RIGHT){
            std::cout << "STATE after_slide=";
            printCurrentStates();
        }

        // Update exit rates
        updateExitRates(selectedSubChainIndex_);
        std::cout << "RATES total_exit_rate_after_update=" << totalExitRate_ << "\n";
        std::cout << "RATES subchain_exit_rates_after_update=";
        for (double rate : subChainExitRates_) {
            std::cout << rate << ",";
        }
        std::cout << "\n";

        // Record event time + event
        pastTimes_.push_back(currentTime_);
        pastEvents_.push_back(PastEvent{selectedSubChainIndex_, nextEdit});
        std::cout << "STEP_RESULT event_recorded_time=" << currentTime_ << "\n";
        ++step;
    }
}

void MultiTFMarkovChain::printCurrentStates() const{
    for (const auto& subChainPtr : currentChain_) {
        const State& state = subChainPtr->getCurrentState();
        int stateID = state.getStateID(topo_.getNumSides(), topo_.getNumProteins());
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

    std::cout << "Done. Wrote CSV to: " << filename << "\n";

}
