//
// Created by jan on 22/10/18.
//

#include "src/LFNS/mpi/MpiTags.h"
#include "src/base/IoUtils.h"
#include "src/particle_filter/ParticleFilter.h"
#include "src/options/LFNSOptions.h"
#include "src/base/Utils.h"
#include "ABCSetup.h"
#include "src/ABC/ABCMpi.h"
#include "src/ABC/ABCWorker.h"


static std::string model_summary_suffix = "model_summary";
static std::stringstream model_summary_stream;

options::LFNSOptions lfns_options;

int my_rank;
int num_tasks;

namespace bmpi = boost::mpi;

double computeDist(const std::vector<double> &theta, ABCSetup abc_setup);

void runMaster(ABCSetup &lfns_setup);

void runWorker(ABCSetup &abc_setup);

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

        ABCSetup lfns_setup(lfns_options, my_rank);
//        _abc_setup = &lfns_setup;
        lfns_setup.setUp();
        if (my_rank == 0) { runMaster(lfns_setup); }
        else { runWorker(lfns_setup); }
    } catch (const std::exception &e) {
        std::cerr << "Failed to run LFNS, exception thrown:\n\t" << e.what() << std::endl;
        return 0;
    }
}


void runMaster(ABCSetup &lfns_setup) {

    abc::ABCMpi lfns(lfns_setup.abc_settings, num_tasks);
    lfns.setSampler(lfns_setup.prior, lfns_setup.density_estimation, lfns_setup.rng);
    lfns.setLogParams(lfns_setup.sampler_settings.getLogParams());
    if (!lfns_setup.abc_settings.previous_log_file.empty()) {
        lfns.resumeRum(lfns_setup.abc_settings.previous_log_file);
    }
    lfns.runABC();
}

void runWorker(ABCSetup &abc_setup) {

    if (my_rank == 1) {
        abc_setup.printSettings(model_summary_stream);
        std::cout << model_summary_stream.str();

        std::string model_summary_file_name = base::IoUtils::appendToFileName(
                abc_setup.abc_settings.output_file,
                model_summary_suffix);
        std::ofstream model_summary_file_stream(model_summary_file_name.c_str());
        model_summary_file_stream << model_summary_stream.str();
        model_summary_file_stream.close();
    }

    abc::LogLikelihodEvalFct_ptr cmp_fct_ptr = std::make_shared<abc::LogLikelihodEvalFct>(
            std::bind(computeDist, std::placeholders::_1, abc_setup));

    abc::ABCWorker worker(my_rank, abc_setup.full_models.front()->getUnfixedParamteters().size(), cmp_fct_ptr);
    worker.setSampler(abc_setup.prior, abc_setup.density_estimation, abc_setup.rng);
    worker.setLogParams(abc_setup.sampler_settings.getLogParams());

    for (int i = 0; i < abc_setup.particle_filters.size(); i++) {
        abc_setup.simulators[i]->addStoppingCriterion(worker.getStoppingFct());
        abc_setup.particle_filters[i]->addStoppingCriterion(worker.getStoppingFct());
        if (abc_setup.particle_filter_settings.use_premature_cancelation) {
            abc_setup.particle_filters[i]->setThresholdPtr(worker.getEpsilonPtr());
        }
    }
    if (abc_setup.particle_filter_settings.use_premature_cancelation) {
        abc_setup.mult_like_eval.setThresholdPtr(worker.getEpsilonPtr());
    }
    worker.run();
}

particle_filter::LogLikelihodEvalFct_ptr getDistFct() {

}

typedef std::vector<double> Times;
typedef std::vector<double> State;

double computeDist(const std::vector<double> &theta, ABCSetup abc_setup) {
    double t = 0.0;
    State latentstate(abc_setup.full_models[0]->dynamics->getNumSpecies(), 0.0);
    State measurement(abc_setup.full_models[0]->measurement_model->getNumMeasurements(), 0.0);

    abc_setup.simulators[0]->reset(latentstate, t);
    abc_setup.full_models[0]->setParameter(theta);

    Trajectory latent_states;
    latent_states.reserve(abc_setup.data_vec[0].size());
    Trajectory measurements;
    measurements.reserve(abc_setup.data_vec[0].size());

    abc_setup.full_models[0]->initial_value_provider->computeInitialState(&latentstate, &t);
    abc_setup.simulators[0]->reset(latentstate, t);

    int num_sim = abc_setup.particle_filter_settings.H;
    double distance = 0;
    for (int i = 0; i < num_sim; i++) {
//        std::cout << "average taken over " << num_sim << " simulations" << std::endl;
        int count = 0;
        for (double sim_time : abc_setup.times_vec[0]) {
            try {
                abc_setup.simulators[0]->simulate(sim_time);

            } catch (mu::ParserError &e) {
                std::ostringstream os;
                os << "Parser error for expression : " << e.GetExpr() << std::endl;
                os << "Message:  " << e.GetMsg() << "\n";
                os << "Formula:  " << e.GetExpr() << "\n";
                os << "Token:    " << e.GetToken() << "\n";
                os << "Position: " << e.GetPos() << "\n";
                os << "Errc:     " << e.GetCode() << "\n";
                throw std::runtime_error(os.str());
            }
            abc_setup.full_models[0]->measurement_model->computeMeasurement(&measurement, latentstate, t);
            latent_states.push_back(latentstate);
            measurements.push_back(measurement);
            TrajectorySet &tmp = abc_setup.data_vec.back();
            Trajectory &tmp_2 = tmp[0];
            std::vector<double> &tmp_3 = tmp_2[count];

            distance += std::pow(measurement[0] - tmp_3[0], 2);
            count++;
        }
    }
    distance = distance / num_sim;

//    std::cout << "returning " << distance << " with " << num_sim << " simulations" << std::endl;
    return distance;

}

