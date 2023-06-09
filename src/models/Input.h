//
// Created by jan on 6/9/23.
//

#ifndef LFNS_INPUT_H
#define LFNS_INPUT_H

#include <string>
#include <memory>

namespace models {

    struct InputData {
        explicit InputData(std::string input_name) : input_name(input_name) {}

        InputData(const InputData &rhs) : input_name(rhs.input_name) {}

        virtual ~InputData()= default;

        InputData &operator=(const InputData &rhs) {
            if (this == &rhs) { return *this; }
            *this = rhs;
            this->input_name = rhs.input_name;
            return *this;
        }

        std::string input_name;
    };
    typedef std::shared_ptr<InputData> InputData_ptr;

    class Input {
    public:
        explicit Input(InputData input_data) : input_name(input_data.input_name),
                                      parameter_index(-1) {};

        virtual void evaluateInput(std::vector<double> &modified_parameter, const double *state, double t) = 0;

        virtual std::vector<double> getDisContTimes() = 0;

        std::string input_name;
        int parameter_index;
    };

    typedef std::shared_ptr<Input> Input_ptr;

} // models

#endif //LFNS_INPUT_H
