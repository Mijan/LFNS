//
// Created by jan on 22/10/18.
//

#include "src/LFNS/mpi/MpiTags.h"
#include "src/base/IoUtils.h"
#include "src/particle_filter/ParticleFilter.h"
#include "src/LFNS/mpi/LFNSWorker.h"
#include "src/options/LFNSOptions.h"
#include "src/base/Utils.h"
#include "LFNSSetup.h"
#include "src/LFNS/mpi/LFNSMpi.h"


static std::string model_summary_suffix = "model_summary";
static std::stringstream model_summary_stream;

LFNSSetup lfns_setup;

options::LFNSOptions lfns_options;

int my_rank;
int num_tasks;

namespace bmpi = boost::mpi;

void runMaster();

void runWorker();

int main(int argc, char *argv[]) {
    bmpi::environment env;
    bmpi::communicator world;

    my_rank = world.rank();
    num_tasks = world.size();

    if (num_tasks == 1) {
        std::cerr << "When using mpi version of LFNS at least two processes must be started!" << std::endl;
        abort();
    }

    try {
        lfns_options.handleCommandLineOptions(argc, argv);

        std::cout << "Config file: " << lfns_options.config_file_name << std::endl;
        lfns_setup.setUp(lfns_options);
        if (my_rank == 0) { runMaster(); }
        else { runWorker(); }
    } catch (const std::exception &e) {
        std::cerr << "Failed to run LFNS, exception thrown:\n\t" << e.what() << std::endl;
        return 0;
    }
}


void runMaster() {

    lfns::mpi::LFNSMpi lfns(lfns_setup.settings, lfns_setup.rng, num_tasks);
    if (!lfns_setup.settings.previous_log_file.empty()) { lfns.resumeRum(lfns_setup.settings.previous_log_file); }
    lfns.runLFNS();
}

void runWorker() {

    if (my_rank == 1) {
        std::size_t max_num_traj = 0;
        for (TrajectorySet &data : lfns_setup.data_vec) {
            max_num_traj = max_num_traj > data.size() ? max_num_traj : data.size();
        }

        lfns_setup.settings.num_used_trajectories = std::min((int) max_num_traj,
                                                             lfns_setup.settings.num_used_trajectories);
        lfns_setup.settings.print(model_summary_stream);
        lfns_setup.full_models.front()->printInfo(model_summary_stream);

        std::cout << model_summary_stream.str();

        std::string model_summary_file_name = base::IoUtils::appendToFileName(lfns_setup.settings.output_file,
                                                                              model_summary_suffix);
        std::ofstream model_summary_file_stream(model_summary_file_name.c_str());
        model_summary_file_stream << model_summary_stream.str();
        model_summary_file_stream.close();
    }

    lfns::mpi::LFNSWorker worker(my_rank, lfns_setup.full_models.front()->getUnfixedParamteters().size(),
                                 lfns_setup.mult_like_eval.getLogLikeFun());

    for (int i = 0; i < lfns_setup.particle_filters.size(); i++) {
        lfns_setup.simulators[i]->addStoppingCriterion(worker.getStoppingFct());
        lfns_setup.particle_filters[i]->addStoppingCriterion(worker.getStoppingFct());
        if (lfns_setup.settings.use_premature_cancelation) {
            lfns_setup.particle_filters[i]->setThresholdPtr(worker.getEpsilonPtr());
        }
    }
    worker.run();
}
