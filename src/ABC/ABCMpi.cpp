//
// Created by jan on 05/10/18.
//

#include <iostream>
#include "ABCMpi.h"
#include "../LFNS/mpi/RequestQueue.h"
#include "../LFNS/mpi/MpiRequest.h"
#include "../sampler/DpGmmSampler.h"
#include "../sampler/EllipsoidSampler.h"
#include "../sampler/GaussianSampler.h"
#include "../sampler/KernelDensityEstimation.h"
#include "../sampler/KernelSupportEstimation.h"
#include "../sampler/RejectionSupportSampler.h"
#include "../sampler/UniformSampler.h"
#include "../base/IoUtils.h"


#include <boost/serialization/export.hpp>

BOOST_SERIALIZATION_ASSUME_ABSTRACT(sampler::Sampler);
BOOST_SERIALIZATION_ASSUME_ABSTRACT(sampler::KernelSampler);

BOOST_CLASS_EXPORT_GUID(sampler::Sampler, "Sampler");
BOOST_CLASS_EXPORT_GUID(sampler::KernelSampler, "KernelSampler");
BOOST_CLASS_EXPORT_GUID(sampler::DensityEstimation, "DensityEstimation");
BOOST_CLASS_EXPORT_GUID(sampler::DpGmmSampler, "DpGmmSampler");
BOOST_CLASS_EXPORT_GUID(sampler::EllipsoidSampler, "EllipsoidSampler");
BOOST_CLASS_EXPORT_GUID(sampler::KernelDensityEstimation, "KernelDensityEstimation");
BOOST_CLASS_EXPORT_GUID(sampler::KernelSupportEstimation, "KernelSupportEstimation");
BOOST_CLASS_EXPORT_GUID(sampler::RejectionSupportSampler, "RejectionSupportSampler");
BOOST_CLASS_EXPORT_GUID(sampler::UniformSampler, "UniformSampler");
BOOST_CLASS_EXPORT_GUID(sampler::GaussianSampler, "GaussianSampler");

namespace abc {

    ABCMpi::ABCMpi(ABCSettings &abc_settings, int num_tasks)
            : ABC(abc_settings), _num_tasks(num_tasks) {}

    ABCMpi::~ABCMpi() {}

    void ABCMpi::runABC() {
        bool lfns_terminate = false;

        int m = 0;

        lfns::mpi::RequestQueue queue;
        _initializeQueue(queue);
        if (!_resume_run) {
            _logger.abcStarted(m, _epsilon);
            _samplePrior(queue);
            _live_points.writeToFile(_settings.output_file, "live_points_0");
        } else {
            m = _logger.iterationNumber();
            _epsilon = _logger.lastEpsilon();
            _logger.lfnsResume(m, _epsilon);
        }
        while (!lfns_terminate) {
            m++;
            _logger.iterationStarted(m);

            _sampleNewPoint(queue);
            lfns_terminate = _postIteration();
        }
        _logger.abcTerminated();
        for (int i = 0; i < _num_tasks; i++) { world.send(i, lfns::mpi::INSTRUCTION, lfns::mpi::DIETAG); }
    }

    void ABCMpi::_samplePrior(lfns::mpi::RequestQueue &queue) {
        while (_live_points.numberParticles() < _settings.N) {
            std::queue<std::size_t> &finished_tasks = queue.getFinishedProcessess();
            while (!finished_tasks.empty()) {
                queue.addRequest(finished_tasks.front(), _num_parameters, true);
                finished_tasks.pop();
            }

            if (queue.firstParticleFinished()) {
                double l = queue.getFirstLikelihood();
                const std::vector<double> &theta = queue.getFirstTheta();
                _logger.thetaSampled(theta, queue.getFirstSamplingClocks());
                _logger.distanceComputed(l);
                _live_points.push_back(theta, l);
                _logger.particleAccepted(theta, l, queue.getFirstParticleClocks(), queue.getFirstUsedProcess());
                queue.clearFirstParticle();
            }
        }
        _logger.logIterationStats();
        queue.stopPendingRequests();
    }

    void ABCMpi::_initializeQueue(lfns::mpi::RequestQueue &queue) {
        for (std::size_t rank = 1; rank < _num_tasks; rank++) { queue.addRequest(rank, _num_parameters, true); }
    }

    void ABCMpi::_sampleNewPoint(lfns::mpi::RequestQueue &queue) {
        _sampler->updateAcceptanceRate(_logger.lastAcceptanceRate());
        for (int j = 0; j < _settings.r; j++) {
            const lfns::LFNSParticle &particle = _live_points.removeHighestPartcile();
        }

        std::cout <<"getting epsilon" << std::endl;
        _epsilon = _live_points.getHighestLogLikelihood();
        std::cout << "new epsilon: " << _epsilon << std::endl;

        _updateEpsilon(_epsilon);
        _logger.epsilonUpdated(_epsilon);

        time_t tic = clock();
        _sampler->updateLiveSamples(_live_points);
        time_t toc = clock();
        _logger.samplerUpdated(*_sampler, toc - tic);
        _updateSampler();

//        _live_points = lfns::LiveParticleSet();
        while (_live_points.numberParticles() < _settings.N) {
            std::queue<std::size_t> &finished_tasks = queue.getFinishedProcessess();

            while (!finished_tasks.empty()) {
                queue.addRequest(finished_tasks.front(), _num_parameters);
                finished_tasks.pop();
            }
            while (queue.firstParticleFinished() && _live_points.numberParticles() < _settings.N) {
                double l = queue.getFirstLikelihood();
                const std::vector<double> &theta = queue.getFirstTheta();
                _logger.thetaSampled(theta, queue.getFirstSamplingClocks());
                _logger.distanceComputed(l);
//                std::cout <<"received " << l << std::endl;
                if (l < _epsilon) {
                    _live_points.push_back(theta, l);
                    _logger.particleAccepted(theta, l, queue.getFirstParticleClocks(), queue.getFirstUsedProcess());
                }
                queue.clearFirstParticle();
            }
        }
        std::cout << "finished... " << std::endl;
        queue.stopPendingRequests();
    }


    void ABCMpi::_updateEpsilon(double epsilon) {
        for (int rank = 1; rank < _num_tasks; rank++) {
            world.send(rank, lfns::mpi::INSTRUCTION, lfns::mpi::UPDATE_EPSILON);
            world.send(rank, lfns::mpi::EPSILON, epsilon);
        }
    }


    void ABCMpi::_updateSampler() {
        for (int rank = 1; rank < _num_tasks; rank++) {
            world.send(rank, lfns::mpi::INSTRUCTION, lfns::mpi::UPDATE_SAMPLER);

            std::stringstream stream;
            _sampler->getSerializedSampler(stream);

            std::size_t sampler_size = stream.str().size();
            world.send(rank, lfns::mpi::SAMPLER_SIZE, sampler_size);

            std::string sampler_string = stream.str();
            world.send(rank, lfns::mpi::SAMPLER, sampler_string.c_str(), sampler_size);
        }
    }
}