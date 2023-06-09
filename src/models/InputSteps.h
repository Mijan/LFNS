//
// Created by jan on 6/9/23.
//

#ifndef LFNS_INPUTSTEPS_H
#define LFNS_INPUTSTEPS_H


#include <vector>
#include <string>

namespace models {

    struct StepsData{
        StepsData(std::vector<double> time_pts, std::vector<double> input_strength) : time_pts(
                time_pts), input_strength(input_strength){}

        StepsData(const StepsData &rhs) : time_pts(
                rhs.time_pts), input_strength(rhs.input_strength) {}

        StepsData &operator=(const StepsData &rhs) {
            if (this == &rhs) { return *this; }
            time_pts = rhs.time_pts;
            input_strength = rhs.input_strength;
            return *this;
        }

        std::vector<double> time_pts;
        std::vector<double> input_strength;


    };
    class InputSteps {
        InputSteps(StepsData input_data);

        virtual ~InputSteps();

        std::vector<double> getDisContTimes();

        bool pulseActive(double t);

        std::vector<double> pulse_beginnings;
        std::vector<double> pulse_ends;
        std::string input_name;
        double _input_strength;
        int parameter_index;
    };

} // models

#endif //LFNS_INPUTSTEPS_H
