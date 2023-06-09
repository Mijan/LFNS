//
// Created by jan on 6/9/23.
//

#ifndef LFNS_INPUTSTEPS_H
#define LFNS_INPUTSTEPS_H


#include <vector>
#include <string>
#include "Input.h"

namespace models {

    struct StepData : public InputData {
        explicit StepData(std::string input_name, std::vector<double> time_pts, std::vector<double> input_strengths)
                : InputData(input_name), time_pts(time_pts), input_strengths(input_strengths) {}

        StepData(const StepData &rhs) : InputData(rhs), time_pts(rhs.time_pts), input_strengths(rhs.input_strengths) {}

        ~StepData(){}

        StepData &operator=(const StepData &rhs) {
            if (this == &rhs) { return *this; }
            *this = rhs;
            this->time_pts = rhs.time_pts;
            this->input_strengths = rhs.input_strengths;
            return *this;
        }

        std::vector<double> time_pts;
        std::vector<double> input_strengths;
    };

    class InputSteps : public Input {
    public:
        explicit InputSteps(StepData step_data);

        virtual ~InputSteps();

        std::vector<double> getDisContTimes() override;

        void evaluateInput(std::vector<double> &modified_parameter, const double *state, double t) override;

    protected:
        std::vector<double> _time_pts;
        std::vector<double> _input_strengths;
    };

} // models

#endif //LFNS_INPUTSTEPS_H
