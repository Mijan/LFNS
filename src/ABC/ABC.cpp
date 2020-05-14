//
// Created by jan on 05/10/18.
//

#include <iostream>
#include "ABC.h"
#include "../LFNS/LiveParticleSet.h"
#include "../LFNS/DeadParticleSet.h"
#include "../LFNS/LFNSParticle.h"

namespace abc {
    ABC::ABC(ABCSettings &abc_settings)
            : _settings(abc_settings), _live_points(), _sampler(nullptr),
              _logger(abc_settings), _resume_run(false), _epsilon(1e10), _epsilon_ptr(&_epsilon),
              _num_parameters(-1) {}

    ABC::~ABC() {}

    void ABC::resumeRum(std::string previous_log_file) {
        try {
            _logger.readFromFile(previous_log_file);
        } catch (const std::runtime_error &e) {
            std::stringstream ss;
            ss << "Failed to read previous population from file " << previous_log_file << ":\n\t";
            ss << e.what() << std::endl;
            throw std::runtime_error(ss.str());
        }
        int it_nbr = _logger.iterationNumber();

        size_t pos = previous_log_file.find("_log_file.txt");
        std::string file_start;
        file_start.assign(previous_log_file.begin(), previous_log_file.begin() + pos);

        std::stringstream ss;
        ss << file_start << "_live_points_" << it_nbr << ".txt";
        std::string live_points_file = ss.str();

        _live_points.readFromFile(live_points_file);


        _resume_run = true;

    }

    double *ABC::getPointerToThreshold() { return _epsilon_ptr; }

    bool ABC::checkIfInitialized(std::ostream &os) {
        if (!_sampler.get()) {
            os << "LFNS_Sampler has not been set!" << std::endl;
            return false;
        }
        return true;
    }

    void
    ABC::setSampler(sampler::Sampler_ptr prior, sampler::DensityEstimation_ptr density_estimation, base::RngPtr rng) {
        _sampler = std::make_shared<lfns::LFNSSampler>(prior, density_estimation, rng);
        _num_parameters = prior->getSamplerDimension();
    }

    void ABC::setLogParams(std::vector<int> log_params) { _sampler->setLogParams(log_params); }

    void ABC::setThresholdPointer(double *epsilon_ptr) { _epsilon_ptr = epsilon_ptr; }

    bool ABC::_postIteration() {
        _logger.logIterationStats();
        _logger.logIterationResults();
        std::stringstream ss;
        ss << "live_points_" << _logger.iterationNumber();
        _live_points.writeToFile(_settings.output_file, ss.str());
        _logger.writeToFile();
        bool lfns_terminate = _testTermination(_logger.lastEpsilon());
        return lfns_terminate;
    }

    bool ABC::_testTermination(double curr_dist) {
        return std::log(curr_dist) < _settings.log_termination;
    }
}