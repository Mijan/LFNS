//
// Created by jan on 28/10/18.
//

#ifndef ABC_ABCSETUP_H
#define ABC_ABCSETUP_H


#include "src/simulator/Simulator.h"
#include "src/particle_filter/ParticleFilter.h"
#include "src/models/FullModel.h"
#include "src/LFNS/LFNSSettings.h"
#include "src/options/LFNSOptions.h"
#include "src/LFNS/LFNS.h"
#include "src/io/IoSettings.h"
#include "src/particle_filter/ParticleFilterSettings.h"
#include "src/sampler/SamplerSettings.h"
#include "src/particle_filter/MultLikelihoodEval.h"
#include "src/io/ConfigFileInterpreter.h"
#include "GeneralSetup.h"
#include "src/ABC/ABCSettings.h"

typedef std::vector<double> Times;
typedef std::vector<std::vector<double>> Trajectory;
typedef std::vector<Trajectory> TrajectorySet;

class ABCSetup : public GeneralSetup {
public:
    ABCSetup(options::LFNSOptions options, int process_nbr = 1);

    virtual ~ABCSetup();

    std::vector<particle_filter::ParticleFilter_ptr> particle_filters;
    particle_filter::MultLikelihoodEval mult_like_eval;

    particle_filter::ParticleFilterSettings particle_filter_settings;
    sampler::Sampler_ptr prior;
    sampler::DensityEstimation_ptr density_estimation;
    sampler::SamplerSettings sampler_settings;
    abc::ABCSettings abc_settings;

    std::vector<Times> times_vec;
    std::vector<TrajectorySet> data_vec;
    double threshold = -DBL_MAX;


    void setUp();

    void printSettings(std::ostream &os) override;

private:
    options::LFNSOptions _lfns_options;

    void _readSettingsfromFile() override;

    std::vector<std::string> _readExperiments() override;

    TrajectorySet
    _createData(int num_outputs, std::string experiment, particle_filter::ParticleFilterSettings &settings);

    Times _createDataTimes(std::string experiment, particle_filter::ParticleFilterSettings &settings);

    particle_filter::ParticleFilterSettings _readParticleFilterSettings();

    sampler::DensityEstimation_ptr _createDensityEstimation(abc::ABCSettings lfns_settings,
                                                            sampler::SamplerSettings settings);

    sampler::Sampler_ptr _createPrior(abc::ABCSettings lfns_settings, sampler::SamplerSettings settings);

    sampler::SamplerSettings _readSamplerSettings();

    abc::ABCSettings _readABCSettings();
};


#endif //ABC_ABCSETUP_H
