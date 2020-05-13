//
// Created by jan on 22/10/18.
//

#include "ABCWorker.h"

namespace abc {
    ABCWorker::ABCWorker(std::size_t my_rank, int num_parameters,
                         abc::LogLikelihodEvalFct_ptr log_likelihood_evaluation) :
            _my_rank(my_rank), _num_parameters(num_parameters), _particle(new double[num_parameters]),
            _epsilon(-DBL_MAX), _sampler_size(1), _stopping_flag_request(new bmpi::request()),
            _mpi_stopping_criterion(std::make_shared<MPIStoppingCriterion>(_stopping_flag_request)),
            _log_likelihood_evaluation(log_likelihood_evaluation), _stop_iteration(false),
            _sampler(nullptr) {}

    ABCWorker::~ABCWorker() {
        delete[] _particle;
        delete _stopping_flag_request;
    }


    StoppingFct_ptr ABCWorker::getStoppingFct() {
        return std::make_shared<StoppingFct>(
                std::bind(&MPIStoppingCriterion::processStopped, _mpi_stopping_criterion.get()));
    }

    void ABCWorker::run() {


        if (!_sampler.get()) {
            std::stringstream ss;
            ss << "Tried to run LFNS without initializing it:\n\t" << "LFNS_Sampler has not been set!" << std::endl;
            throw std::runtime_error(ss.str());
        }

        lfns::mpi::MPI_INSTRUCTION instruction(lfns::mpi::INSTRUCTION);

        bool process_terminated = false;
        while (!process_terminated) {
            world.recv(0, lfns::mpi::INSTRUCTION, instruction);
            switch (instruction) {

                case lfns::mpi::DIETAG:
                    std::cout << "This is process " << _my_rank << " , I am exiting." << std::endl;
                    process_terminated = true;
                    break;

                case lfns::mpi::LIKELIHOOD_RECOMPU:
                    _computeLikelihood();
                    break;

                case lfns::mpi::SAMPLE_CONSTR_PRIOR:
                    _sampleConstrPrior();
                    break;

                case lfns::mpi::SAMPLE_PRIOR:
                    _samplePrior();
                    break;

                case lfns::mpi::UPDATE_EPSILON:
                    _updateEpsilon();
                    break;

                case lfns::mpi::UPDATE_SAMPLER:
                    _updateSampler();
                    break;

                case lfns::mpi::PREPARE_STOPPING:
                    _prepareStoppingFlag();
                    break;

                default:
                    break;
            }
        }

        return;
    }


    double *ABCWorker::getEpsilonPtr() { return &_epsilon; }

    void ABCWorker::_computeLikelihood() {
        world.recv(0, lfns::mpi::PARTICLE, _particle, _num_parameters);
        std::vector<double> parameter(_particle, _particle + _num_parameters);
        double log_likelihood = (*_log_likelihood_evaluation)(parameter);
        world.send(0, lfns::mpi::LIKELIHOOD_RECOMPU, log_likelihood);
    }


    void ABCWorker::_sampleConstrPrior() {
        time_t tic = clock();
        const std::vector<double> &parameter = _sampler->sampleConstrPrior();
        time_t toc = clock();
        time_t sampling_time = toc - tic;
        double log_likelihood = (*_log_likelihood_evaluation)(parameter);
//        if (log_likelihood < _epsilon) {
            world.send(0, lfns::mpi::INSTRUCTION, lfns::mpi::PARTICLE_ACCEPTED);
            world.send(0, lfns::mpi::LIKELIHOOD_RECOMPU, log_likelihood);
            world.send(0, lfns::mpi::PARTICLE, parameter.data(), _num_parameters);
            world.send(0, lfns::mpi::CLOCKS_SAMPLING, sampling_time);
//        } else {
//            world.send(0, lfns::mpi::INSTRUCTION, lfns::mpi::PARTICLE_REJECTED);
//        }
    }

    void ABCWorker::_samplePrior() {
        time_t tic = clock();
        std::vector<double> parameter = _sampler->samplePrior();
        time_t toc = clock();
        double log_likelihood = (*_log_likelihood_evaluation)(parameter);
        world.send(0, lfns::mpi::INSTRUCTION, lfns::mpi::PARTICLE_ACCEPTED);
        world.send(0, lfns::mpi::LIKELIHOOD_RECOMPU, log_likelihood);
        world.send(0, lfns::mpi::PARTICLE, parameter.data(), _num_parameters);;
        world.send(0, lfns::mpi::CLOCKS_SAMPLING, toc - tic);
    }

    void ABCWorker::_updateSampler() {
        world.recv(0, lfns::mpi::SAMPLER_SIZE, _sampler_size);
        char sampler_char_ptr[_sampler_size];
        world.recv(0, lfns::mpi::SAMPLER, sampler_char_ptr, _sampler_size);
        std::string sampler_str(sampler_char_ptr, sampler_char_ptr + _sampler_size);
        std::stringstream stream(sampler_str);
        _sampler->updateSerializedSampler(stream);
        _sampler->getDensityEstimation()->updateLogLikelihoodFct(_log_likelihood_evaluation);
    }

    void ABCWorker::_prepareStoppingFlag() {
        _stop_iteration = false;
        *_stopping_flag_request = world.irecv(0, lfns::mpi::STOP_SIMULATION, _stop_iteration);
        _mpi_stopping_criterion->updateRequest(_stopping_flag_request);
    }


    void
    ABCWorker::setSampler(sampler::Sampler_ptr prior, sampler::DensityEstimation_ptr density_estimation,
                          base::RngPtr rng) {
        _sampler = std::make_shared<lfns::LFNSSampler>(prior, density_estimation, rng);
    }

    void ABCWorker::setLogParams(std::vector<int> log_params) { _sampler->setLogParams(log_params); }


    void ABCWorker::_updateEpsilon() { world.recv(0, lfns::mpi::EPSILON, _epsilon); }
}

