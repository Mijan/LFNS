//
// Created by jan on 05/10/18.
//

#include <iostream>
#include <iomanip>
#include "ABCLogger.h"
#include "../base/Utils.h"
#include "../base/MathUtils.h"
#include "../base/IoUtils.h"
#include "../LFNS/LFNSSampler.h"

namespace abc {

    ABCLogger::ABCLogger(ABCSettings settings) : _remaining_required_particles_iteration(settings.N),
                                                 _num_samples_iteration(0), _num_samples_particle(0),
                                                 _num_accepted_iteration(0),
                                                 _print_interval(settings.print_interval),
                                                 _acceptance_info_print_interval(
                                                         settings.acceptance_info_print_interval), _particle_tic(0),
                                                 _algorithm_tic(0), _iteration_tic(0),
                                                 _output_file_name(settings.output_file), _acceptance_rates(),
                                                 _epsilons(), _iteration_nbrs(), _seconds_for_iteration(),
                                                 _sampling_seconds_for_iteration() {}

    ABCLogger::~ABCLogger() {}

    void ABCLogger::abcStarted(int m, double epsilon) {
        std::cout << "\n\nStarting ABC algorithm at initial iteration " << m << " and initial epsilon " << epsilon
                  << std::endl << std::endl;
        _algorithm_tic = clock();
    }

    void ABCLogger::lfnsResume(int m, double epsilon) {
        std::cout << "Resume ABC algorithm after initial iteration " << m << " and last epsilon " << epsilon
                  << std::endl << std::endl;
    }

    void ABCLogger::iterationStarted(int m) {
        std::cout << "\n\nStarting iteration " << m << std::endl;
        _iteration_nbrs.push_back(m);
        _num_accepted_iteration = 0;
        _num_samples_iteration = 0;
        _iteration_tic = clock();
    }

    void ABCLogger::deadPointAdded(const lfns::LFNSParticle &particle) {
        _remaining_required_particles_iteration++;
    }

    void ABCLogger::epsilonUpdated(double epsilon) {
        std::cout << "New epsilon: " << epsilon << std::endl;
        _epsilons.push_back(epsilon);
    }

    void ABCLogger::samplerUpdated(lfns::LFNSSampler &sampler, time_t clocks_for_sampler_update) {
        std::cout << "Sampler updated: ";
        sampler.writeToStream(std::cout);
        std::cout << std::endl;
    }

    void ABCLogger::thetaSampled(const std::vector<double> &theta, time_t clocks_sampling) {
        _num_samples_particle++;
        _num_samples_iteration++;
        _particle_tic = clock();
    }

    void ABCLogger::distanceComputed(double likelihood) {}

    void ABCLogger::particleAccepted(const std::vector<double> &theta, double l) {
        particleAccepted(theta, l, clock() - _particle_tic);
    }


    void ABCLogger::particleAccepted(const std::vector<double> &theta, double l, time_t clocks_particle,
                                     int origin_process) {

        _num_accepted_iteration++;
        _remaining_required_particles_iteration--;

        if (_num_accepted_iteration % _print_interval == 0) {
            std::cout << "i = " << _num_accepted_iteration << ", log likelihood = " << l
                      << ", acceptance probability = " << 1.0 / _num_samples_particle << ", time for particle = "
                      << clocks_particle / (double) CLOCKS_PER_SEC;
            if (origin_process > 0) { std::cout << ", process: " << origin_process; }
            std::cout << std::endl;
        }

        if (_num_accepted_iteration % _acceptance_info_print_interval == 0 ||
            _remaining_required_particles_iteration == 0) { _printAcceptanceInfo(); }

        _num_samples_particle = 0;
    }

    void ABCLogger::logIterationResults() {
        double acceptance_rate = _num_accepted_iteration / (double) _num_samples_iteration;

        clock_t toc2 = clock();
        std::cout << std::endl << "\nIteration " << _iteration_nbrs.back()
                  << " needed " << base::Utils::getTimeFromClocks(toc2 - _iteration_tic) << std::endl;
        std::cout << "Total LF-NS time so far: " << base::Utils::getTimeFromClocks(toc2 - _algorithm_tic) << std::endl;
        std::cout << "Acceptance rate: " << acceptance_rate << std::endl;
    }

    void ABCLogger::logIterationStats() {
        double acceptance_rate = _num_accepted_iteration / (double) _num_samples_iteration;
        _acceptance_rates.push_back(acceptance_rate);

        clock_t toc2 = clock();

        double sec = ((double) (toc2 - _iteration_tic) / CLOCKS_PER_SEC);
        _seconds_for_iteration.push_back(sec);

    }

    void ABCLogger::abcTerminated() {
        time_t toc = clock();
        std::cout << "\n\nABC algorithm successfully terminated!" << std::endl;
        std::cout << "Total ABC time: " << base::Utils::getTimeFromClocks(toc - _algorithm_tic) << std::endl
                  << std::endl;
    }

    void ABCLogger::_printAcceptanceInfo() {
        double acceptance_rate = ((double) _num_accepted_iteration)
                                 / ((double) _num_samples_iteration);
        std::cout << std::endl << std::endl << "Particle acceptance rate: "
                  << acceptance_rate << std::endl;

        clock_t toc2 = clock();
        std::cout << "Time needed so far for this iteration: " << base::Utils::getTimeFromClocks(toc2 - _iteration_tic)
                  <<
                  std::endl;

        double eta = ((toc2 - _iteration_tic) / ((double) _num_accepted_iteration))
                     * _remaining_required_particles_iteration;
        std::cout << "Estimated time remaining: " << base::Utils::getTimeFromClocks(eta) <<
                  std::endl << std::endl;
    }

    int ABCLogger::iterationNumber() {
        std::cout << "returning the last of " << _iteration_nbrs.size() << std::endl;
        return _iteration_nbrs.back();
    }


    double ABCLogger::lastEpsilon() { return _epsilons.back(); }

    double ABCLogger::lastAcceptanceRate() { return _acceptance_rates.back(); }

    void ABCLogger::writeToFile() {
        std::string log_file_name =
                base::IoUtils::appendToFileName(_output_file_name,
                                                "log_file");
        std::ofstream log_file_file(log_file_name.c_str());
        if (!log_file_file.is_open()) {
            std::cerr << "error opening file "
                      << log_file_name.c_str()
                      << " for writing run log.. run log could not be saved!!"
                      << std::endl << std::endl;
            return;
        }
        log_file_file << std::setw(3) << "i" << "\t" << std::setw(15) << std::setprecision(8) << "epsilon"
                      << "\t" << std::setw(20) << std::setprecision(8) << "Acceptance rate" << "\t"
                      << std::setw(12) << std::setprecision(8) << "seconds" << "\t" << std::setw(17) << std::endl;

        for (size_t i = 0; i < _iteration_nbrs.size(); i++) {
            int num_simulation = _iteration_nbrs[i];
            double epsilon = _epsilons[i];
            double acceptance_rate = _acceptance_rates[i];
            double seconds = _seconds_for_iteration[i];


            log_file_file << std::setw(3) << num_simulation << "\t" << std::setw(15) << std::setprecision(6) << epsilon
                          << "\t" << std::setw(20) << std::setprecision(6) << acceptance_rate << "\t"
                          << std::setw(12) << std::setprecision(8) << seconds << std::endl;
        }
        log_file_file.close();
        std::cout << "Log wrote into " << log_file_name.c_str() << std::endl;
    }

    void ABCLogger::readFromFile(std::string previous_log_file_name) {
        std::string log_file_name = previous_log_file_name;
        std::ifstream log_file(log_file_name.c_str());
        if (!log_file.is_open()) {
            std::stringstream ss;
            ss << "error opening file " << log_file_name.c_str() << " for reading logs of previous run!"
               << std::endl;
            throw std::runtime_error(ss.str());
        }

        std::string line;
        getline(log_file, line); // read the title
        std::istringstream iss(line);


        int sim_nbr;
        double epsilon;
        double acceptance_rate;
        double seconds;

        while (!log_file.eof()) {

            std::getline(log_file, line);
            if (!line.empty()) {
                iss.clear();
                iss.str(line);


                std::string num_simulation_str;
                iss >> num_simulation_str;
                sim_nbr = std::stoi(num_simulation_str.c_str());
                _iteration_nbrs.push_back(sim_nbr);

                std::string epsilon_str;
                iss >> epsilon_str;
                epsilon = std::stod(epsilon_str.c_str());
                _epsilons.push_back(epsilon);

                std::string acceptance_rate_string;
                iss >> acceptance_rate_string;
                acceptance_rate = std::stod(acceptance_rate_string.c_str());
                _acceptance_rates.push_back(acceptance_rate);
                iss >> seconds;
                _seconds_for_iteration.push_back(seconds);
            }
        }

        log_file.close();
        std::cout << "Log read from " << log_file_name.c_str() << std::endl;
    }
}
