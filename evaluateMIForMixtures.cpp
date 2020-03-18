//
// Created by jan on 2/18/20.
//

#include <string>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <DPGMMEstimator.h>
#include "src/base/IoUtils.h"
#include "src/base//EigenMatrices.h"
#include "src/sampler/Sampler.h"
#include "src/LFNS/LFNSSettings.h"
#include "src/LFNS/LFNS.h"
#include "src/LFNS/seq/LFNSSeq.h"
#include "src/sampler/DpGmmSampler.h"
#include "src/sampler/RejectionSupportSampler.h"
#include "src/base/MultivariateNormal.h"
#include "src/sampler/UniformSampler.h"

namespace fs = boost::filesystem;
namespace po = boost::program_options;

std::string config_file_name;
std::string output_file_name;
po::variables_map vm;

int N = 200;

sampler::DensityEstimation_ptr sampler_XY;
sampler::Sampler_ptr sampler_X;
sampler::Sampler_ptr sampler_Y;

bool H_computation = false;

DP_GMM::GaussMixtureComponentSet readMixture(std::string file_name);

DP_GMM::GaussMixtureComponentSet getXMixtures(DP_GMM::GaussMixtureComponentSet &mixtures);

DP_GMM::GaussMixtureComponentSet getYMixtures(DP_GMM::GaussMixtureComponentSet &mixtures);

double logLike(const std::vector<double> &sample);

class MixutreSampler : public sampler::DensityEstimation {
public:
    MixutreSampler(base::RngPtr rng, DP_GMM::GaussMixtureComponentSet components, sampler::SamplerData data)
            : sampler::DensityEstimation(rng, data), comps(components), _dist(0, 1) {}

    virtual ~MixutreSampler() {};
    DP_GMM::GaussMixtureComponentSet comps;
    base::UniformRealDistribution _dist;


    virtual std::vector<double> &sample() {
        const DP_GMM::GaussMixtureComponentPtr comp;
        DP_GMM::GaussMixtureComponentSet::iterator it = comps.begin();
        double u = _dist(*_rng);
        double w = 0.0;
        w += it->get()->comp_weight;
        while (w < u) {
            it++;
            w += it->get()->comp_weight;
        }

        const base::EiVector &mean = it->get()->mean;
        const base::EiMatrix &decomp_cov = it->get()->decomposed_cov;

        base::MultivariateNormal::mvnormRndWithDecomposedVar(_rng, _sample.size(), mean, decomp_cov, &_sample);
        return _sample;
    };

    virtual double getLogLikelihood(const std::vector<double> &sample) {
        double likelihood = 0.0;

        DP_GMM::GaussMixtureComponentSet::const_iterator it;
        for (it = comps.begin(); it != comps.end(); it++) {
            const DP_GMM::GaussMixtureComponentPtr comp = *it;
            likelihood += comp->comp_weight
                          * base::MultivariateNormal::mvnormPrepared(
                    sample.size(), sample, comp->mean, comp->precision,
                    1.0 / comp->precision_det);
        }
        return log(likelihood);
    }

    virtual void updateTransformedDensitySamples(const base::EiMatrix &transformed_samples) {}

    virtual void sampleTransformed(base::EiVector &trans_sample) {}

    virtual double getTransformedLogLikelihood(const base::EiVector &trans_sample) {}
};


int main(int argc, char **argv) {
    try {
        po::options_description desc;
        desc.add_options()("help", "produce help message")("config_file",
                                                           po::value<std::string>(&config_file_name),
                                                           "Config file. A config file must always be provided!")(
                "output_file,O", po::value<std::string>(&output_file_name), "Output file")(
                "entory,h", po::value<bool>(&H_computation), "Compute entropy")(
                "num_runs,N", po::value<int>(&N), "number of LFNS particles");

        po::positional_options_description p;
        p.add("config_file", -1);

        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 0;
        }


        DP_GMM::GaussMixtureComponentSet mixtures = readMixture(config_file_name);
        DP_GMM::GaussMixtureComponentSet::iterator it = mixtures.begin();
        int num_dim = it->get()->mean.size();

        base::RngPtr rng = std::make_shared<base::RandomNumberGenerator>(1 * time(NULL));


        sampler::SamplerData data(num_dim);
        sampler_XY = std::make_shared<MixutreSampler>(rng, mixtures, data);


        if (!H_computation) {
            DP_GMM::GaussMixtureComponentSet mixtures_X = getXMixtures(mixtures);
            DP_GMM::GaussMixtureComponentSet mixtures_Y = getYMixtures(mixtures);

            sampler::SamplerData dataX(num_dim - 1);
            sampler_X = std::make_shared<MixutreSampler>(rng, mixtures_X, dataX);
            sampler::SamplerData dataY(num_dim - 1);
            sampler_Y = std::make_shared<MixutreSampler>(rng, mixtures_Y, dataY);
        }

        lfns::LFNSSettings settings;
        settings.log_termination = -3.9;
        settings.uniform_prior = true;
        settings.output_file = output_file_name;
        settings.r = N / 2;
        settings.N = N;

        lfns::LogLikelihodEvalFct_ptr log_likelihood_evaluation = std::make_shared<lfns::LogLikelihodEvalFct>(&logLike);
        lfns::seq::LFNSSeq lfns_seq(settings, log_likelihood_evaluation);


        sampler::DpGmmSamplerData dpgmm_data(num_dim);
        dpgmm_data.num_dp_iterations = 50;
        sampler::DpGmmSampler_ptr dpgmm_sampler = std::make_shared<sampler::DpGmmSampler>(rng, dpgmm_data);


        sampler::RejectionSamplerData rej_data(data);
        rej_data.rejection_quantile = 0.1;

        sampler::DensityEstimation_ptr density_estimation_ptr = std::make_shared<sampler::RejectionSupportSampler>(rng,
                                                                                                                   dpgmm_sampler,
                                                                                                                   rej_data);

        sampler::UniformSamplerData uniform_data(num_dim);
        for (int i = 0; i < num_dim; i++) {
            std::pair<double, double> bound = {0, 3};;
            uniform_data.bounds[i] = bound;
        }
//        sampler::Sampler_ptr prior = std::make_shared<sampler::UniformSampler>(rng, uniform_data);

        lfns_seq.setSampler(sampler_XY, sampler_XY, rng);

        lfns_seq.runLFNS();

    } catch (const std::exception &e) {
        std::cerr << "Failed to run MI estimation, exception thrown:\n\t" << e.what() << std::endl;
        return 0;
    }
}


DP_GMM::GaussMixtureComponentSet readMixture(std::string file_name) {

    DP_GMM::GaussMixtureComponentSet mixture_componentes;
    std::string means_file_name = base::IoUtils::appendToFileName(file_name.c_str(), "means");
    std::string cov_file_name = base::IoUtils::appendToFileName(file_name.c_str(), "covariances");
    std::string weights_file_name = base::IoUtils::appendToFileName(file_name.c_str(), "weights");

    std::vector<double> means_vector = base::IoUtils::readColVector(means_file_name);
    std::vector<double> cov_vector = base::IoUtils::readColVector(cov_file_name);
    std::vector<double> weights_vector = base::IoUtils::readColVector(weights_file_name);

    std::size_t num_components = weights_vector.size();
    std::size_t num_dimensions = means_vector.size() / num_components;

    int mean_index = 0;
    int cov_index = 0;
    for (int i = 0; i < num_components; i++) {
        double weight = weights_vector[i];

        base::EiVector mean(num_dimensions);
        base::EiMatrix cov_(num_dimensions, num_dimensions);
        for (int j = 0; j < num_dimensions; j++) {
            mean(j) = means_vector[mean_index++];
            for (int k = 0; k < num_dimensions; k++) {
                cov_(j, k) = cov_vector[cov_index++];
            }
        }


        base::EiMatrix precision = cov_.inverse();
        double prec_det = precision.determinant();
        DP_GMM::GaussMixtureComponentPtr comp = std::make_shared<DP_GMM::GaussMixtureComponent>(mean, cov_, precision,
                                                                                                prec_det, weight);

        mixture_componentes.insert(comp);


    }
    return mixture_componentes;
}


DP_GMM::GaussMixtureComponentSet getXMixtures(DP_GMM::GaussMixtureComponentSet &mixtures) {
    DP_GMM::GaussMixtureComponentSet mixture_componentes;


    DP_GMM::GaussMixtureComponentSet::iterator it;

    for (it = mixtures.begin(); it != mixtures.end(); it++) {
        DP_GMM::GaussMixtureComponentPtr orig_comp = *it;

        base::EiVector new_mean(1);
        new_mean(0) = orig_comp->mean(0);

        base::EiMatrix new_cov(1, 1);
        new_cov(0, 0) = orig_comp->cov(0, 0);

        base::EiMatrix new_prec(1, 1);
        new_prec = new_cov.inverse();

        double det = new_prec.determinant();

        DP_GMM::GaussMixtureComponentPtr new_comp = std::make_shared<DP_GMM::GaussMixtureComponent>(new_mean, new_cov,
                                                                                                    new_prec, det,
                                                                                                    orig_comp->comp_weight);
        mixture_componentes.insert(new_comp);
    }
    return mixture_componentes;
}

DP_GMM::GaussMixtureComponentSet getYMixtures(DP_GMM::GaussMixtureComponentSet &mixtures) {
    DP_GMM::GaussMixtureComponentSet mixture_componentes;


    DP_GMM::GaussMixtureComponentSet::iterator it;

    for (it = mixtures.begin(); it != mixtures.end(); it++) {
        DP_GMM::GaussMixtureComponentPtr orig_comp = *it;
        int num_dim = orig_comp->mean.size() - 1;
        base::EiVector new_mean(num_dim);
        base::EiMatrix new_cov(num_dim, num_dim);

        for (int j = 0; j < num_dim; j++) {
            new_mean(j) = orig_comp->mean(j + 1);
            for (int k = 0; k < num_dim; k++) {
                new_cov(j, k) = orig_comp->cov(j + 1, k + 1);
            }
        }


        base::EiMatrix new_prec(new_mean.size(), new_mean.size());
        new_prec = new_cov.inverse();

        double det = new_prec.determinant();

        DP_GMM::GaussMixtureComponentPtr new_comp = std::make_shared<DP_GMM::GaussMixtureComponent>(new_mean, new_cov,
                                                                                                    new_prec, det,
                                                                                                    orig_comp->comp_weight);
        mixture_componentes.insert(new_comp);
    }
    return mixture_componentes;
}


double logLike(const std::vector<double> &sample) {
    double log_pxy = sampler_XY->getLogLikelihood(sample);

    double result = 0.0;
    if (H_computation) {
        result = log(-log_pxy);
    } else {
        std::vector<double> sample_X = {sample[0]};
        std::vector<double> sample_Y(sample.begin() + 1, sample.end());

        double log_px = sampler_X->getLogLikelihood(sample_X);
        double log_py = sampler_Y->getLogLikelihood(sample_Y);

        result = log(log_pxy - log_px - log_py);
    }

    if (std::isnan(result)) { result = -DBL_MAX; }
    return result;

}