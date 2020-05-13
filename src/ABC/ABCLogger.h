//
// Created by jan on 05/10/18.
//

#ifndef ABC_ABCLOGGER_H
#define ABC_ABCLOGGER_H

#include <vector>
#include "../LFNS/LFNSParticle.h"
#include "../LFNS/LiveParticleSet.h"
#include "../sampler/DensityEstimation.h"
#include "ABCSettings.h"
#include "../LFNS/LFNSSampler.h"

namespace abc {
    class ABCLogger {

    public:
        ABCLogger(ABCSettings settings);

        virtual ~ABCLogger();

        void abcStarted(int i, double d);

        void lfnsResume(int m, double epsilon);

        void iterationStarted(int i);

        void deadPointAdded(const lfns::LFNSParticle &particle);

        void epsilonUpdated(double d);

        void samplerUpdated(lfns::LFNSSampler &sampler, time_t clocks_for_sampler_update);

        void thetaSampled(const std::vector<double> &vector, time_t clocks_sampling);

        void distanceComputed(double likelihood);

        void particleAccepted(const std::vector<double> &theta, double l);

        void
        particleAccepted(const std::vector<double> &theta, double l, time_t clocks_particle, int origin_process = 0);

        void logIterationResults();

        void logIterationStats();

        void abcTerminated();

        void writeToFile();

        void readFromFile(std::string previous_log_file_name);

        int iterationNumber();

        double lastEpsilon();

        double lastAcceptanceRate();

    private:
        int _remaining_required_particles_iteration;

        int _num_samples_iteration;
        int _num_samples_particle;
        int _num_accepted_iteration;

        int _print_interval;
        int _acceptance_info_print_interval;

        time_t _particle_tic;
        time_t _algorithm_tic;
        time_t _iteration_tic;

        std::string _output_file_name;
        std::vector<double> _acceptance_rates;
        std::vector<double> _epsilons;
        std::vector<int> _iteration_nbrs;
        std::vector<double> _seconds_for_iteration;
        std::vector<double> _sampling_seconds_for_iteration;

        void _printAcceptanceInfo();
    };
}


#endif //ABC_ABCLOGGER_H
