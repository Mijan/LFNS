//
// Created by jan on 6/9/23.
//

#include "InputSteps.h"
#include "../base/MathUtils.h"

namespace models {
    InputSteps::InputSteps(StepData step_data)
            : Input(step_data), _time_pts(step_data.time_pts), _input_strengths(step_data.input_strengths) {}

    InputSteps::~InputSteps(){};

    std::vector<double> InputSteps::getDisContTimes(){
        return _time_pts;
    }

    void InputSteps::evaluateInput(std::vector<double> &modified_parameter, const double *state, double t){
        if (t < _time_pts.front() || t > _time_pts.back()) { return; }
        int pulse_index = base::MathUtils::findIndexSmallerThan(_time_pts, t);
        double input_strength = _input_strengths[pulse_index];
        modified_parameter[parameter_index] += input_strength;
    }

} // models