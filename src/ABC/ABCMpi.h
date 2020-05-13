//
// Created by jan on 05/10/18.
//

#ifndef ABC_ABCMPI_H
#define ABC_ABCMPI_H

#include "ABC.h"
#include "ABCSettings.h"
#include "../LFNS/mpi/MpiTags.h"
#include "../LFNS/mpi/RequestQueue.h"


namespace abc {

        class ABCMpi : public ABC {

        public:
            ABCMpi(ABCSettings &abc_settings, int num_tasks);

            virtual ~ABCMpi();

            void runABC() override;

        private:

            void _updateSampler(sampler::DensityEstimation_ptr density_estimation);

            void _requestDistance(std::vector<double> &theta);


            bmpi::communicator world;
            int _num_tasks;

            void _samplePrior(lfns::mpi::RequestQueue &queue);

            void _sampleNewPoint(lfns::mpi::RequestQueue &queue);

            void _updateEpsilon(double epsilon);

            void _updateSampler();

            void _initializeQueue(lfns::mpi::RequestQueue &queue);
        };
}


#endif //ABC_ABCMPI_H
