//
// Created by jan on 25/02/19.
//

#include "MCMCSampler.h"
#include "../sampler/UniformSampler.h"
#include "../sampler/GaussianSampler.h"

namespace mcmc {

    MCMCSampler::MCMCSampler(MCMCSettings &mcmc_settings, sampler::SamplerSettings &settings, base::RngPtr rng) : _rng(
            rng), _prior(), _uniform_prior(true), _log_params() {

        std::vector<std::string> unfixed_params = settings.param_names;
        sampler::SamplerData sampler_data(unfixed_params.size());
        sampler_data.bounds = settings.getBounds(unfixed_params);
        _log_params = settings.getLogParams(unfixed_params);

        for (int i : _log_params) {
            sampler_data.bounds[i].first = std::log10(sampler_data.bounds[i].first);
            sampler_data.bounds[i].second = std::log10(sampler_data.bounds[i].second);
        }

        if (mcmc_settings.uniform_prior) {
            sampler::UniformSamplerData uni_data(sampler_data);
            _prior = std::make_shared<sampler::UniformSampler>(_rng, uni_data);
        }

        switch (mcmc_settings.kernel_type) {
            case GAUSS  : {
                sampler::NormalSamplerData normal_data(sampler_data);
//                normal_data.cov = base::EiMatrix::Identity(sampler_data.size(), sampler_data.size()) * 0.5;
                // LacGfp empirical cov
//                normal_data.cov = base::EiMatrix(18, 18);
//                normal_data.cov << 2.8984,0.034621,-0.11615,0.26603,0.30667,-0.040102,0.0070074,-0.019816,0.045401,0.01132,-0.083132,0.0066597,0.27667,0.023078,0.0057327,0.019157,-0.00072689,-0.041338,0.034621,3.1332,0.27378,0.1444,0.20314,0.25602,-0.26536,-0.022124,0.0413,-0.11854,-0.11194,0.043168,0.052433,0.18077,0.01999,0.009493,-0.0060854,-0.05951,-0.11615,0.27378,2.6753,0.16836,0.00056151,-0.12722,0.053945,-0.1354,-0.095675,-0.20838,-0.030453,0.055123,0.065637,0.044414,0.0047575,0.0093226,0.0056295,-0.06448,0.26603,0.1444,0.16836,3.0143,-0.20209,-0.03289,-0.41164,-0.021276,-0.010929,-0.0383,0.029486,0.027214,0.094801,0.0014995,0.020142,0.043267,0.003163,-0.088019,0.30667,0.20314,0.00056151,-0.20209,3.0311,0.047984,0.1796,0.10605,0.02792,-0.050132,-0.13376,0.0012422,-0.20902,-0.2137,0.043533,0.025743,-0.018142,-0.016881,-0.040102,0.25602,-0.12722,-0.03289,0.047984,2.3576,0.034261,-0.23213,0.027338,-0.0060026,-0.020074,0.0038118,0.0018008,0.082008,0.0078792,-0.011494,-0.0054384,0.0080668,0.0070074,-0.26536,0.053945,-0.41164,0.1796,0.034261,2.5235,0.14018,0.036635,-0.045962,-0.029365,-0.015269,0.076415,-0.17283,0.0093818,-0.0045368,-0.014027,0.059184,-0.019816,-0.022124,-0.1354,-0.021276,0.10605,-0.23213,0.14018,0.28977,0.032699,-0.0039291,0.020712,0.020706,0.027948,-0.071643,-0.023731,0.032048,-0.00085015,-0.0011675,0.045401,0.0413,-0.095675,-0.010929,0.02792,0.027338,0.036635,0.032699,0.49468,0.52604,-0.31041,-0.01455,0.30301,0.027846,0.020705,0.0034548,-0.0093992,0.021154,0.01132,-0.11854,-0.20838,-0.0383,-0.050132,-0.0060026,-0.045962,-0.0039291,0.52604,1.2743,0.084369,-0.019819,0.40193,0.11174,-0.006992,0.0079898,-0.0049224,0.004432,-0.083132,-0.11194,-0.030453,0.029486,-0.13376,-0.020074,-0.029365,0.020712,-0.31041,0.084369,0.73801,0.022751,-0.17046,-0.019566,-0.028497,-0.0072201,0.013251,-0.015828,0.0066597,0.043168,0.055123,0.027214,0.0012422,0.0038118,-0.015269,0.020706,-0.01455,-0.019819,0.022751,0.12867,0.06264,0.081322,0.010697,-0.034236,0.003579,-0.080874,0.27667,0.052433,0.065637,0.094801,-0.20902,0.0018008,0.076415,0.027948,0.30301,0.40193,-0.17046,0.06264,2.1431,0.13051,0.031175,0.0091722,0.0017192,-0.051599,0.023078,0.18077,0.044414,0.0014995,-0.2137,0.082008,-0.17283,-0.071643,0.027846,0.11174,-0.019566,0.081322,0.13051,1.2621,0.027061,-0.024678,0.0037068,-0.055334,0.0057327,0.01999,0.0047575,0.020142,0.043533,0.0078792,0.0093818,-0.023731,0.020705,-0.006992,-0.028497,0.010697,0.031175,0.027061,0.056259,0.0040635,-0.016931,-0.015743,0.019157,0.009493,0.0093226,0.043267,0.025743,-0.011494,-0.0045368,0.032048,0.0034548,0.0079898,-0.0072201,-0.034236,0.0091722,-0.024678,0.0040635,0.1158,0.0030803,-0.082624,-0.00072689,-0.0060854,0.0056295,0.003163,-0.018142,-0.0054384,-0.014027,-0.00085015,-0.0093992,-0.0049224,0.013251,0.003579,0.0017192,0.0037068,-0.016931,0.0030803,0.01399,-0.008925,-0.041338,-0.05951,-0.06448,-0.088019,-0.016881,0.0080668,0.059184,-0.0011675,0.021154,0.004432,-0.015828,-0.080874,-0.051599,-0.055334,-0.015743,-0.082624,-0.008925,0.20816;
                // lotka-voltera empirical cov
                normal_data.cov = base::EiMatrix(3, 3);
                normal_data.cov
                        << 0.000355993778619239, 0.000194657465316225, 0.000207652176430899, 0.000194657465316225, 0.000288027332982342, 0.000177855727729757, 0.000207652176430899, 0.000177855727729757, 0.000399517173562744;
                sampler::GaussianSampler_ptr gauss_kernel = std::make_shared<sampler::GaussianSampler>(_rng,
                                                                                                       normal_data);

                _kernel_sampler = gauss_kernel;
                break;
            }
        }
    }


    const std::vector<double> &MCMCSampler::samplePrior() { return _scaleUpSample(_prior->sample()); }

    const std::vector<double> &MCMCSampler::sampleKernel(const std::vector<double> &kernel_center) {
        std::vector<double> scale_center = kernel_center;
        scale_center = _scaleDownSample(scale_center);
        std::vector<double> &sample = _kernel_sampler->sample(scale_center);
        while (!_prior->isSampleFeasible(sample)) {
            sample = _kernel_sampler->sample(scale_center);
        }
        return _scaleUpSample(sample);
    }


    std::vector<double> &MCMCSampler::_scaleUpSample(std::vector<double> &sample) {
        if (_log_params.empty()) { return sample; }
        else {
            for (int &index : _log_params) {
                sample[index] = std::pow(10, sample[index]);
            }
            return sample;
        }
    }

    std::vector<double> &MCMCSampler::_scaleDownSample(std::vector<double> &sample) {
        if (_log_params.empty()) { return sample; }
        else {
            for (int &index : _log_params) {
                sample[index] = std::log10(sample[index]);
            }
            return sample;
        }
    }

}