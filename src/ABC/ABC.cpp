//
// Created by jan on 05/10/18.
//

#include <iostream>
#include "ABC.h"
#include "../LFNS/LiveParticleSet.h"
#include "../LFNS/DeadParticleSet.h"
#include "../LFNS/LFNSParticle.h"
#include "../base/IoUtils.h"
#include "../base/MathUtils.h"

namespace abc {
    ABC::ABC(ABCSettings &abc_settings)
            : _settings(abc_settings), _live_points(),        _sampler(nullptr),
              _logger(abc_settings), _resume_run(false), _epsilon(-DBL_MAX), _epsilon_ptr(&_epsilon),
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

        std::string dead_points_file = file_start + "_dead_points.txt";

        std::stringstream ss;
        ss << file_start << "_live_points_" << it_nbr << ".txt";
        std::string live_points_file = ss.str();

        if (base::IoUtils::doesFileExists(dead_points_file) && base::IoUtils::doesFileExists(live_points_file)) {
            _live_points.readFromFile(live_points_file);
        } else {
            std::string posterior_file_name = file_start + "_posterior.txt";
            std::string posterior_log_like_file_name = file_start + "_posterior_log_likelihoods.txt";

            if (base::IoUtils::doesFileExists(posterior_file_name) &&
                base::IoUtils::doesFileExists(posterior_log_like_file_name)) {
                _live_points.readFromFile(posterior_file_name);
                if (_live_points.numberParticles() == it_nbr * _settings.r + _settings.N) {
                    for (int i = 0; i < it_nbr * _settings.r; i++) {
                        const lfns::LFNSParticle &particle = _live_points.removeLowestPartcile();
                    }
                } else {
                    std::cerr << "Previous posterior files found, but number of particles does not match!";
                    std::cerr << " log file indicates that the last iteration was iteration " << it_nbr << " and N="
                              << _settings.N << ", and r = " << _settings.r;
                    std::cerr << ". Thus posterior is expected to have " << it_nbr * _settings.r << " + " << _settings.N
                              << "=" << it_nbr * _settings.r + _settings.N << " particles, but posterior has "
                              << _live_points.numberParticles() << " particles." << std::endl;
                    _live_points = lfns::LiveParticleSet();
                }
            } else {
                std::cerr << "Neither dead points and live points files " << dead_points_file << " and "
                          << live_points_file << " or posterior file " << posterior_file_name << " provided."
                          << std::endl;
            }
        }

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
        std::cout <<"writing to file..." << std::endl;
        _live_points.writeToFile(_settings.output_file, ss.str());
        std::cout << "live points written " << std::endl;
        _logger.writeToFile();
        std::cout <<"written to file" << std::endl;
        bool lfns_terminate = _testTermination(_logger.lastEpsilon());
        std::cout << "and termination is " << lfns_terminate << std::endl;
        return lfns_terminate;
    }

    bool ABC::_testTermination(double curr_dist) {
        std::cout <<"testing " << curr_dist << " smaller than " << _settings.log_termination << std::endl;
        return std::log(curr_dist) < _settings.log_termination;
    }
}