//
// Created by jan on 05/10/18.
//

#ifndef ABC_ABC_H
#define ABC_ABC_H


#include <functional>
#include "ABCSettings.h"
#include "../LFNS/LiveParticleSet.h"
#include "ABCLogger.h"
#include "../sampler/DensityEstimation.h"
#include "../LFNS/LFNSSampler.h"
#include <iostream>

namespace abc {
    // TODO find better place for typedef
    typedef std::function<double(const std::vector<double> &)> LogLikelihodEvalFct;
    typedef std::shared_ptr<LogLikelihodEvalFct> LogLikelihodEvalFct_ptr;

    class ABC {

    public:
        explicit ABC(ABCSettings &abc_settings);

        virtual ~ABC();

        virtual void runABC() = 0;

        void resumeRum(std::string previous_log_file);

        double *getPointerToThreshold();

        bool checkIfInitialized(std::ostream &os);

        void
        setSampler(sampler::Sampler_ptr prior, sampler::DensityEstimation_ptr density_estimation, base::RngPtr rng);

        void setLogParams(std::vector<int> log_params);

        void setThresholdPointer(double * epsilon_ptr);

    protected:
        ABCSettings _settings;

        lfns::LiveParticleSet _live_points;

        lfns::LFNSSampler_ptr _sampler;

        ABCLogger _logger;
        bool _resume_run;
        double *_epsilon_ptr;
        double _epsilon;

        int _num_parameters;

        bool _testTermination(double curr_dist);

        bool _postIteration();
    };
};


#endif //ABC_ABC_H
