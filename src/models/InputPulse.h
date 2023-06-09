//
// Created by jan on 26/10/18.
//

#ifndef LFNS_INPUTPULSE_H
#define LFNS_INPUTPULSE_H

#include <vector>
#include <string>
#include <set>
#include "../base/Utils.h"
#include "Input.h"

namespace models {
    struct PulseData : public InputData {
        PulseData(double pulse_period_, double pulse_strenght_, double pulse_duration_, int num_pulses_,
                  std::string pulse_inpt_name_, double starting_time_) : InputData(pulse_inpt_name_), pulse_period(
                pulse_period_), pulse_strenght(pulse_strenght_), pulse_duration(pulse_duration_),
                                                                         num_pulses(num_pulses_),
                                                                         starting_time(starting_time_) {}


        PulseData(const PulseData &rhs) :
                InputData(rhs), pulse_period(rhs.pulse_period), pulse_strenght(rhs.pulse_strenght),
                pulse_duration(rhs.pulse_duration),
                num_pulses(rhs.num_pulses),
                starting_time(rhs.starting_time) {}

        PulseData &operator=(const PulseData &rhs) {
            if (this == &rhs) { return *this; }
            *this = rhs;
            pulse_period = rhs.pulse_period;
            pulse_strenght = rhs.pulse_strenght;
            pulse_duration = rhs.pulse_duration;
            num_pulses = rhs.num_pulses;
            starting_time = rhs.starting_time;
            return *this;
        }

        double pulse_period;
        double pulse_strenght;
        double pulse_duration;
        int num_pulses;
        double starting_time;
    };

    class InputPulse : public Input {
    public:
        explicit InputPulse(PulseData input_data);

        virtual ~InputPulse();

        std::vector<double> getDisContTimes() override;

        void evaluateInput(std::vector<double> &modified_parameter, const double *state, double t) override;

        std::vector<double> pulse_beginnings;
        std::vector<double> pulse_ends;
        double _input_strength;

    protected:
        bool _pulseActive(double t);

    };

    struct InputPulses {

        InputPulses() : input_pram_indices(), modified_parameter(), pulses() {}

        virtual ~InputPulses() {}

        std::vector<double> getDisContTime() {
            std::vector<double> times;
            for (InputPulse &pulse: pulses) { base::Utils::addOnlyNew(times, pulse.getDisContTimes()); }
            std::sort(times.begin(), times.end());
            return times;
        }

        std::vector<double> modified_parameter;
        std::set<int> input_pram_indices;
        std::vector<InputPulse> pulses;
    };
}


#endif //LFNS_INPUTPULSE_H
