//
// Created by jan on 6/9/23.
//

#ifndef LFNS_INPUT_H
#define LFNS_INPUT_H

#include <string>
#include <memory>

namespace models {

    class Input {
    public:
        Input(std::string input_name, int parameter_index) : input_name(input_name),
                                                             parameter_index(parameter_index) {};

        virtual void evaluateInput(std::vector<double> &modified_parameter, const double *state, double t) = 0;

        virtual std::vector<double> getDisContTimes() = 0;

        std::string input_name;
        int parameter_index;
    };
    typedef std::shared_ptr<Input> Input_ptr;

} // models

#endif //LFNS_INPUT_H
