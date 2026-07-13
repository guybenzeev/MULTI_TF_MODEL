// MultiTFMarkovChain.cpp
#include "MultiTFMarkovChain.h"
#include "FenwickTree.h"
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <chrono>


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

int MultiTFMarkovChain::selectSubChainFromTree(double& threshold) {
    int index = ratesTree_.findIndex(threshold);

    if (index < 0 || index >= length_) {
        return -1;
    }

    double previousCumulativeRate = 0.0;

    if (index > 0) {
        previousCumulativeRate = ratesTree_.prefixSum(index - 1);
    }

    threshold -= previousCumulativeRate;

    return index;
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

        ratesTree_.add(i, rate);

        if(rate > largestExitRate_){
            largestExitRate_ = rate;
        }
    }
}

void MultiTFMarkovChain::updateExitRates(int editHeadIndex, int proteinWidth) {
    
    const int largestWidth = topo_.getLargestProteinWidth();
    const int largestRadius = topo_.getLargestRadius();

    int bottom = editHeadIndex - largestWidth - largestRadius;
    int top = editHeadIndex + proteinWidth + largestRadius;

    if (bottom < 0) {
        bottom = 0;
    }

    if (top > length_) {
        top = length_;
    }

    for(int i = bottom; i < top; ++i){
        double oldRate = subChainExitRates_[i];
        double newRate = currentChain_[i]->getTotalExitRate(currentChain_);
        subChainExitRates_[i] = newRate;
        totalExitRate_ += (newRate - oldRate);
        if(newRate > largestExitRate_){
            largestExitRate_ = newRate;
        }

        ratesTree_.add(i, newRate - oldRate);
    }   

}

void MultiTFMarkovChain::applyEdit(Edit& nextEdit) {
        Edit slidToEdit;
        slidToEdit.type = EditType::BIND_NS;
        int slidingProtein = currentChain_[selectedSubChainIndex_]->getCurrentState().protein;
        /*
        std::cout << "STATE selected_subchain_before_apply="
                  << currentChain_[selectedSubChainIndex_]->getCurrentState().getStateID(topo_.getNumSides(), topo_.getNumProteins())
                  << ", sliding_protein=" << slidingProtein
                  << "\n";
        */
        int width = topo_.getProteinWidth(nextEdit.protein - 1);

        for(int i = 0; i < width; ++i) {
            currentChain_[selectedSubChainIndex_ + i]->applyEdit(nextEdit, selectedSubChainIndex_);
        }

        //CHECK THIS PART FOR SLIDING UPDATES
        if(nextEdit.type == EditType::SLIDE_LEFT) {
            slidToEdit.strandSide = nextEdit.strandSide;
            slidToEdit.protein = slidingProtein;
            for(int i = 0; i < width; ++i) {
                currentChain_[selectedSubChainIndex_ - 1 + i]->applyEdit(slidToEdit, selectedSubChainIndex_-1);
            }
            /*
            std::cout << "SLIDE destination_index=" << (selectedSubChainIndex_ - 1)
                      << ", destination_state="
                      << currentChain_[selectedSubChainIndex_ - 1]->getCurrentState().getStateID(topo_.getNumSides(), topo_.getNumProteins())
                      << "\n";

            */
            updateExitRates(selectedSubChainIndex_ - 1, width);
        }
        else if(nextEdit.type == EditType::SLIDE_RIGHT) {
            slidToEdit.strandSide = nextEdit.strandSide;
            slidToEdit.protein = slidingProtein;
            for(int i = 0; i < width; ++i) {
                currentChain_[selectedSubChainIndex_ + 1 + i]->applyEdit(slidToEdit, selectedSubChainIndex_+1);
            }
            /*
            std::cout << "SLIDE destination_index=" << (selectedSubChainIndex_ + 1)
                      << ", destination_state="
                      << currentChain_[selectedSubChainIndex_ + 1]->getCurrentState().getStateID(topo_.getNumSides(), topo_.getNumProteins())
                      << "\n";
            */
            updateExitRates(selectedSubChainIndex_ + 1, width);
        }
        if(nextEdit.type == EditType::SLIDE_LEFT || nextEdit.type == EditType::SLIDE_RIGHT){
            //std::cout << "STATE after_slide=";
            //printCurrentStates();
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
    /*
    std::cout << "largest protein width: " << topo_.getLargestProteinWidth() << "\n";
    std::ofstream coutLog("data/multitf_cout.log");
    if (!coutLog) {
        throw std::runtime_error("Failed to open cout log file: data/multitf_cout.log");
    }
    CoutRedirect coutRedirect(coutLog);
    */

    pastEvents_.clear();
    pastTimes_.clear();

    largestExitRate_ = 0.0;
    currentTime_ = 0;

    ratesTree_ = FenwickTree(length_);

    initializeExitRates();

    int step = 0;


    //Clear nextEdit?
    while (currentTime_ < runTime) {
        /*std::cout << "\nSTEP " << step << "\n";
        std::cout << "TIME before=" << currentTime_ << "\n";
        std::cout << "RATES total_exit_rate=" << totalExitRate_ << "\n";
        std::cout << "RATES subchain_exit_rates=";
        for (double rate : subChainExitRates_) {
            std::cout << rate << ",";
        }
        std::cout << "\n";
        std::cout << "STATE before=";
        printCurrentStates();
        */
        using Clock = std::chrono::steady_clock; 
        using Microseconds = std::chrono::duration<double, std::micro>; 
        auto loopStart = Clock::now();
        // --------------------------------------------------------- 
        // 1. Sample holding time 
        // ---------------------------------------------------------
        double holdingTime = findHoldingTime(); // sample holding time
        //std::cout << "TIME holding_time=" << holdingTime << "\n";
      
        auto holdingTimeStart = Clock::now();
        currentTime_ += holdingTime;
        //std::cout << "TIME after_holding=" << currentTime_ << "\n";

        auto holdingTimeEnd = Clock::now(); 
        double holdingTimeDuration = Microseconds(holdingTimeEnd - holdingTimeStart).count();
        if (currentTime_ > runTime) {
            //std::cout << "STEP_RESULT exceeded_run_time runTime=" << runTime << "\n";
            break;
        }

        // --------------------------------------------------------- 
        // 2. Generate selection threshold
        // ---------------------------------------------------------

        auto thresholdStart = Clock::now();
        // Select which sub-chain will fire next
        std::uniform_real_distribution<double> dist(0.0, totalExitRate_);
        double threshold = dist(rng_);
        double rawThreshold = threshold;
        auto thresholdEnd = Clock::now();
        double thresholdDuration = Microseconds(thresholdEnd - thresholdStart).count();
        // --------------------------------------------------------- 
        // 3. Select sub-chain
        // ---------------------------------------------------------
   
        //std::cout << "THRESHOLD raw=" << rawThreshold << "\n";
        auto selectSubChainStart = Clock::now();
        selectedSubChainIndex_ = selectSubChainFromTree(threshold);

        //std::cout << "THRESHOLD local_after_subchain_selection=" << threshold << "\n";
        auto selectSubChainEnd = Clock::now();
        double selectSubChainDuration = Microseconds( selectSubChainEnd - selectSubChainStart ).count();

        if (selectedSubChainIndex_ < 0) {
            throw std::runtime_error("selectSubChain failed (threshold out of range?)");
        }
        //std::cout << "SELECTED subchain_index=" << selectedSubChainIndex_
        //          << ", subchain_exit_rate=" << subChainExitRates_[selectedSubChainIndex_]
        //          << "\n";

        // Find and apply the next edit

        // --------------------------------------------------------- 
        // 4. Find next edit
        // ---------------------------------------------------------
        auto findEditStart = Clock::now();
        Edit nextEdit = currentChain_[selectedSubChainIndex_]->findNextEdit(threshold, currentChain_);
        auto findEditEnd = Clock::now();
        double findEditDuration = Microseconds(findEditEnd - findEditStart).count();

        //std::cout << "EDIT selected=" << nextEdit.toString()
        //          << ", protein=" << nextEdit.protein
        //          << "\n";
        if(nextEdit.type == EditType::SLIDE_LEFT || nextEdit.type == EditType::SLIDE_RIGHT){
            //std::cout << "STATE before_slide=";
            //printCurrentStates();
        }

        // --------------------------------------------------------- 
        // 5. Apply edit
        // ---------------------------------------------------------
        auto applyEditStart = Clock::now();
        applyEdit(nextEdit);
        auto applyEditEnd = Clock::now();
        double applyEditDuration = Microseconds(applyEditEnd - applyEditStart).count();
        // Update exit rates

        // --------------------------------------------------------- 
        // 6. Update exit rates
        // ---------------------------------------------------------
        auto updateRatesStart = Clock::now();
        updateExitRates(selectedSubChainIndex_, topo_.getProteinWidth(nextEdit.protein - 1));

        auto updateRatesEnd = Clock::now();
        double updateRatesDuration = Microseconds(updateRatesEnd - updateRatesStart).count();
       // std::cout << "RATES total_exit_rate_after_update=" << totalExitRate_ << "\n";
       // std::cout << "RATES subchain_exit_rates_after_update=";
        //for (double rate : subChainExitRates_) {
        //    std::cout << rate << ",";
        //}
        //std::cout << "\n";

        // Record event time + event
        // --------------------------------------------------------- 
        // 7. Record event 
        // ---------------------------------------------------------
        auto recordEventStart = Clock::now();

        pastTimes_.push_back(currentTime_);
        pastEvents_.push_back(PastEvent{selectedSubChainIndex_, nextEdit});

        auto recordEventEnd = Clock::now();
        double recordEventDuration = Microseconds(recordEventEnd - recordEventStart).count();

        //std::cout << "STEP_RESULT event_recorded_time=" << currentTime_ << "\n";

        // --------------------------------------------------------- 
        // Total loop time 
        // --------------------------------------------------------- 
        auto loopEnd = Clock::now(); 
        double totalLoopDuration = Microseconds(loopEnd - loopStart).count(); 
        /*
        std::cout << "STEP " << step << " timing:\n" << " findHoldingTime: " << 
        holdingTimeDuration << " us\n" << " generateThreshold: " << thresholdDuration 
        << " us\n" << " selectSubChain: " << selectSubChainDuration << " us\n" << " findNextEdit: " 
        << findEditDuration << " us\n" << " applyEdit: " << applyEditDuration << " us\n" << " updateExitRates: " 
        << updateRatesDuration << " us\n" << " recordEvent: " << recordEventDuration << " us\n" << " TOTAL: "
         << totalLoopDuration << " us\n\n";
         */
        ++step;
    }
}

void MultiTFMarkovChain::printCurrentStates() const{
    for (const auto& subChainPtr : currentChain_) {
        const State& state = subChainPtr->getCurrentState();
        int stateID = state.getStateID(topo_.getNumSides(), topo_.getNumProteins());

        if(state.isOccupied()) {
            stateID *= -1; // Mark occupied states with a negative ID
        }

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
   // out << "time,index,edit_type,protein,side\n";
   out << "time,index,edit,protein,side,length=" << length_
    << ",num_sides=" << topo_.getNumSides() << ",num_proteins=" << topo_.getNumProteins() << ",runTime=" << runTime << "\n";

    for (size_t i = 0; i < pastTimes_.size(); ++i) {
        const auto& ev = pastEvents_[i];

        out << std::fixed << std::setprecision(10) << pastTimes_[i] << ","
            << ev.index << ","
            //<< static_cast<int>(ev.edit.type) << ","
            << ev.edit.typeAsString() << ","
            << ev.edit.protein << ","
            << ev.edit.strandSide << "\n";
    }

    out.close();

    std::cout << "Done. Wrote CSV to: " << filename << "\n";

}
