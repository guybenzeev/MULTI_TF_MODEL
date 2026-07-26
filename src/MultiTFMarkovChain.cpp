// MultiTFMarkovChain.cpp
#include "MultiTFMarkovChain.h"
#include "FenwickTree.h"
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <iomanip>

double MultiTFMarkovChain::findHoldingTime() {
    if (totalExitRate_ <= 0.0) {
        throw std::runtime_error("Total exit rate is zero or negative, cannot sample holding time.");
    }

    std::exponential_distribution<double> dist(totalExitRate_);
    return dist(rng_);
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

void MultiTFMarkovChain::initializeExitRates(){
    subChainExitRates_.clear();
    subChainExitRates_.reserve(length_);
    totalExitRate_ = 0.0;
    for(int i = 0; i < length_; ++i){
        double rate = currentChain_[i]->getTotalExitRate(currentChain_);
        subChainExitRates_.push_back(rate);
        totalExitRate_ += rate;

        ratesTree_.add(i, rate);

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
                    i
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

    ratesTree_ = FenwickTree(length_);

    initializeExitRates();

    while (currentTime_ < runTime) {
        // --------------------------------------------------------- 
        // 1. Sample holding time 
        // ---------------------------------------------------------
        double holdingTime = findHoldingTime();
        currentTime_ += holdingTime;

        if (currentTime_ > runTime) {
            break;
        }

        // --------------------------------------------------------- 
        // 2. Generate selection threshold
        // ---------------------------------------------------------

        std::uniform_real_distribution<double> dist(0.0, totalExitRate_);
        double threshold = dist(rng_);

        // --------------------------------------------------------- 
        // 3. Select sub-chain
        // ---------------------------------------------------------
        selectedSubChainIndex_ = selectSubChainFromTree(threshold);

        if (selectedSubChainIndex_ < 0) {
            throw std::runtime_error("selectSubChain failed (threshold out of range?)");
        }

        // --------------------------------------------------------- 
        // 4. Find next edit
        // ---------------------------------------------------------
        Edit nextEdit = currentChain_[selectedSubChainIndex_]->findNextEdit(threshold, currentChain_);

        // --------------------------------------------------------- 
        // 5. Apply edit
        // ---------------------------------------------------------
        applyEdit(nextEdit);

        // --------------------------------------------------------- 
        // 6. Update exit rates
        // ---------------------------------------------------------
        updateExitRates(selectedSubChainIndex_, topo_.getProteinWidth(nextEdit.protein - 1));

        // --------------------------------------------------------- 
        // 7. Record event 
        // ---------------------------------------------------------
        pastTimes_.push_back(currentTime_);
        pastEvents_.push_back(PastEvent{selectedSubChainIndex_, nextEdit});
    }
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
