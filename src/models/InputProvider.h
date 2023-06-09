//
// Created by jan on 6/9/23.
//

#ifndef LFNS_INPUTPROVIDER_H
#define LFNS_INPUTPROVIDER_H

#include <vector>
#include <set>
#include "Input.h"
namespace models {
    struct InputProvider {

        std::vector<double> modified_parameter;
        std::set<int> input_pram_indices;
        std::vector<Input_ptr> inputs;

        std::vector<double> getDisContTime() {
            std::vector<double> times;
            for (Input_ptr input: inputs) { base::Utils::addOnlyNew(times, input->getDisContTimes()); }
            std::sort(times.begin(), times.end());
            return times;
        }
    };
}


#endif //LFNS_INPUTPROVIDER_H
