//
// Created by jan on 05/10/18.
//

#include <cmath>
#include <cfloat>
#include "PosteriorEstimator.h"
#include "../base/MathUtils.h"
#include "../base/IoUtils.h"

namespace lfns {


    PosteriorEstimator::PosteriorEstimator(int N, int r) : _N(N), _r(r), _E_t(), _E_t2(), _e_normalized(),
                                                           _e_log_scale(), _epsilons_normalized(),
                                                           _epsilons_log_scale() {}

    PosteriorQuantitites
    PosteriorEstimator::estimatePosteriorQuantities(LiveParticleSet &live_points, DeadParticleSet &dead_points) {

        _E_t = std::vector<double>(_r + 1, 1.0);
        _E_t2 = std::vector<double>(_r + 1, 1.0);
        for (int j = 1; j <= _r; j++) {
            _E_t[j] = E_t(j);
            _E_t2[j] = E_t2(j);
        }

        _computeE(dead_points);

        PosteriorQuantitites post_quant;
        int m = dead_points.size() / _r;

        std::vector<double> log_delta = getDeltas(m);
        post_quant.log_delta = log_delta;
        double log_final_vol = m * std::log(_E_t[_r]);
        post_quant.log_final_vol = log_final_vol;

        std::vector<double> log_zds = std::vector<double>(m * _r, 0.0);
        double max_log = -DBL_MAX;
        for (int i = 0; i < m * _r; i++) {
            log_zds[i] = dead_points[i].getLogLikelihood() + log_delta[i];
            max_log = max_log > log_zds[i] ? max_log : log_zds[i];
        }
        double log_zd = base::MathUtils::sumOfLog(log_zds, max_log);

        double log_zl = live_points.getLogAverageL() + log_final_vol;
        post_quant.log_zd = log_zd;
        post_quant.log_zl = log_zl;

        post_quant.log_ztot = base::MathUtils::sumOfLog(log_zd, log_zl);

        double log_zd_var = _computeLogZdVar_3(dead_points, log_zd);
        post_quant.log_zd_var = log_zd_var;

        double log_zl_var = _computeLogZlVar(live_points, m);
        post_quant.log_zl_var = log_zl_var;

        double log_cov_term;
        int sign_cov;
        _computeCovTerm(dead_points, live_points, &log_cov_term, &sign_cov);
        double log_ztot_var = base::MathUtils::sumOfLog({log_zl_var, log_zd_var, log_cov_term}, {1, 1, sign_cov});

        post_quant.log_ztot_var = log_ztot_var;

        post_quant.log_live_average = live_points.getLogAverageL();
        post_quant.log_live_variance = live_points.getLogVariance();
        post_quant.log_Ex2 = m * std::log(_E_t2[_r]);

        double log_1 = 0.5 * post_quant.log_ztot_var;
        double log_2 =
                live_points.getLogVariance() - std::log(live_points.numberParticles()) + 2 * m * std::log(_E_t[_r]);
        double log_3 = 0.5 * base::MathUtils::diffOfLog(post_quant.log_ztot_var, log_2);
        post_quant.log_max_std_improvement = base::MathUtils::diffOfLog(log_1, log_3) - post_quant.log_ztot;


        post_quant.log_max_live = live_points.getHighestLogLikelihood();


        double log_average_l = live_points.getLogAverageL();
        double log_var_xm = std::log(_varX(m, _r));

        double v_min = base::MathUtils::sumOfLog({log_var_xm + 2 * log_average_l, log_zd_var, log_cov_term},
                                                 {1, 1, sign_cov});
        post_quant.log_min_var = v_min;
        return post_quant;

    }

    std::vector<double> PosteriorEstimator::getDeltas(int m) {
        std::vector<double> log_x = std::vector<double>(_r * m, 0.0);

        double log_tr = log(E_t(_r));

        for (int j = 0; j < _r; j++) {
            log_x[j] = log(E_t(j + 1));
            for (int k = 1; k < m; k++) {
                log_x[k * _r + j] = log_x[j] + k * log_tr;
            }
        }

        std::vector<double> log_delta = std::vector<double>(_r * m, 0.0);
        log_delta[0] = log(1 - exp(log_x[0]));
        for (int j = 1; j < log_delta.size(); j++) {
            log_delta[j] = log(exp(log_x[j - 1]) - exp(log_x[j]));
        }
        return log_delta;
    }

    double PosteriorEstimator::E_t(int j) { return (_N - j + 1) / (double) (_N + 1); }

    double PosteriorEstimator::E_t2(int j) { return (_N - j + 1) * (_N - j + 2) / (double) ((_N + 1) * (_N + 2)); }




    void PosteriorEstimator::_computeE(DeadParticleSet &dead_points) {

        int m = dead_points.size() / _r;

        _e_normalized.resize(m);
        _e_log_scale.resize(m);

        _epsilons_normalized.resize(m);
        _epsilons_log_scale.resize(m);

        for (int i = 0; i < m; i++) {
            std::vector<double> f(_r + 1, 0.0);
            std::vector<double> epsilon(_r, 0.0);
            double max_log = dead_points[i * _r + (_r - 1)].getLogLikelihood();
            _e_log_scale[i] = max_log;
            _epsilons_log_scale[i] = max_log;
            for (int j = 0; j <= _r; j++) {
                if (j < _r) { epsilon[j] = std::exp(dead_points[i * _r + j].getLogLikelihood() - max_log); }
                if (j == 0) { f[j] = std::exp(dead_points[i * _r].getLogLikelihood() - max_log); }
                if (j > 0 && j < _r) {
                    f[j] = std::exp(dead_points[i * _r + j].getLogLikelihood() - max_log) -
                           std::exp(dead_points[i * _r + j - 1].getLogLikelihood() - max_log);
                }
                if (j == _r) { f[j] = -std::exp(dead_points[i * _r + (_r - 1)].getLogLikelihood() - max_log); }

            }
            _e_normalized[i] = f;
            _epsilons_normalized[i] = epsilon;
        }
    }

    double PosteriorEstimator::_computeLogZlVar(LiveParticleSet &live_points, int m) {
        double log_xm = m * std::log(_E_t[_r]);
        double log_average_l = live_points.getLogAverageL();

        double log_var_l = live_points.getLogVariance(log_average_l) - std::log(live_points.numberParticles());
        double log_var_xm = std::log(_varX(m, _r));

        double sum_2 = log_var_xm + 2 * log_average_l;
        double sum_3 = 2 * log_xm + log_var_l;
        return base::MathUtils::sumOfLog(sum_2, sum_3);

    }

    double PosteriorEstimator::_computeLogZtotVar(double log_zd, double log_zd_var, double log_zl, double log_zl_var,
                                                  LiveParticleSet &live_points, int m) {

        double sum_1 = log_zd_var;
        double sum_2 = log_zl_var;

        double e_tmp = 0;
        std::vector<double> log_facs = std::vector<double>(m, 0.0);
        for (int i = 1; i <= m; i++) {
            double log_fac_2_1 = (i - 1) * log(_E_t2[_r]);

            double log_fac_2_2 = 0;
            for (int j = 0; j < _r; j++) {
                log_fac_2_2 += _e_normalized[i - 1][j] * _E_t[j] * std::pow(_E_t[_r], m - i + 1);
            }
            log_fac_2_2 += _e_normalized[i - 1][_r] * _E_t2[_r] * std::pow(_E_t[_r], m - i);
            log_fac_2_2 = log(log_fac_2_2) + _e_log_scale[i - 1];
            log_facs[i - 1] = log_fac_2_1 + log_fac_2_2;
        }
        e_tmp = base::MathUtils::sumOfLog(log_facs);

        double sum_3_1 = live_points.getLogAverageL() + e_tmp;

        double sum_3_2 = log_zl + log_zd;
        double max_log = std::max<double>(sum_3_1, sum_3_2);

        int sign_3 = exp(sum_3_1 - max_log) - exp(sum_3_2 - max_log) < 0 ? -1 : 1;
        double sum_3 = log(2) + (log(abs(exp(sum_3_1 - max_log) - exp(sum_3_2 - max_log))) + max_log);

        std::vector<double> sums = {sum_1, sum_2, sum_3};
        double max_log_2 = *(std::max_element(sums.begin(), sums.end()));
        double log_ztot_var =
                std::log(exp(sum_1 - max_log_2) + exp(sum_2 - max_log_2) + sign_3 * exp(sum_3 - max_log_2)) + max_log_2;
        return log_ztot_var;
    }

    void PosteriorEstimator::writeToFile(std::string output_file_name, LiveParticleSet &live_points,
                                         DeadParticleSet &dead_points, PosteriorQuantitites &post_quant) {

        int m = dead_points.size() / _N;

        std::string posterior_file_name = base::IoUtils::appendToFileName(output_file_name, "posterior");

        std::ofstream posterior_file(posterior_file_name.c_str());
        if (!posterior_file.is_open()) {
            std::stringstream os;
            os << "error opening file " << posterior_file_name.c_str()
               << " for writing posterior. Posterior could not be saved!!"
               << std::endl;
            throw std::runtime_error(os.str());
        }

        for (std::vector<LFNSParticle>::iterator it = dead_points.begin(); it != dead_points.end(); it++) {
            const std::vector<double> &theta = it->getTheta();
            for (int j = 0; j < theta.size(); j++) {
                double part = theta[j];
                posterior_file << part << " ";
            }
            posterior_file << std::endl;
        }

        for (std::set<LFNSParticle>::iterator it = live_points.begin(); it != live_points.end(); it++) {
            const std::vector<double> &theta = it->getTheta();
            for (int j = 0; j < theta.size(); j++) {
                double part = theta[j];
                posterior_file << part << " ";
            }
            posterior_file << std::endl;
        }
        posterior_file.close();

        std::string likelihood_file_name = base::IoUtils::appendToFileName(output_file_name,
                                                                           "posterior_log_likelihoods");
        std::ofstream likelihood_file(likelihood_file_name.c_str());
        if (!likelihood_file.is_open()) {
            std::stringstream os;
            os << "error opening file " << likelihood_file_name.c_str()
               << " for writing posterior log likelihoods. Posterior log likelihoods could not be saved!!"
               << std::endl;
            throw std::runtime_error(os.str());
        }

        for (std::vector<LFNSParticle>::iterator it = dead_points.begin(); it != dead_points.end(); it++) {
            double log_like = it->getLogLikelihood();
            likelihood_file << log_like << std::endl;
        }

        for (std::set<LFNSParticle>::iterator it = live_points.begin(); it != live_points.end(); it++) {
            double log_like = it->getLogLikelihood();
            likelihood_file << log_like << std::endl;
        }
        likelihood_file.close();

        std::string weight_file_name = base::IoUtils::appendToFileName(output_file_name, "posterior_log_weights");
        std::ofstream weight_file(weight_file_name.c_str());
        if (!weight_file.is_open()) {
            std::stringstream os;
            os << "error opening file " << weight_file_name.c_str()
               << " for writing posterior weights. Posterior weights could not be saved!!"
               << std::endl;
            throw std::runtime_error(os.str());
        }

        std::vector<double> log_deltas = post_quant.log_delta;
        int i = 0;
        for (std::vector<LFNSParticle>::iterator it = dead_points.begin(); it != dead_points.end(); it++) {
            double log_weight = it->getLogLikelihood() + log_deltas[i++];
            weight_file << log_weight << std::endl;
        }


        double log_weight_factor_live_points = post_quant.log_final_vol - std::log(_N);
        for (std::set<LFNSParticle>::iterator it = live_points.begin(); it != live_points.end(); it++) {
            double log_weight = it->getLogLikelihood() + log_weight_factor_live_points;
            weight_file << log_weight << std::endl;
        }
        weight_file.close();


        std::cout << "Posterior written into " << posterior_file_name << std::endl;
    }


    double PosteriorEstimator::_varX(int i, int j) {

        if (j < _r) {
            return std::pow(_E_t2[_r], i - 1) * _E_t2[j] - std::pow(_E_t[_r], 2 * (i - 1)) * std::pow(_E_t[j], 2);
        }
        return std::pow(_E_t2[_r], i) - std::pow(_E_t[_r], 2 * i);
    }

    double PosteriorEstimator::_computeLogZdVar_3(DeadParticleSet &dead_points, double zd) {
        int m_max = dead_points.size() / _r;

        std::vector<double> log_contributions(m_max, 0);
        std::vector<int> signs(m_max, 1);


        for (int m = 1; m <= m_max; m++) {
            std::vector<double> &epsilons = _epsilons_normalized[m - 1];
            double log_part_1 = -DBL_MAX;
            int sign_1 = 1;

            if (m > 1) {

                std::vector<double> log_parts(m - 1, 0.0);
                std::vector<int> signs_parts(m - 1, 1);


                double E_m = 0;
                for (int j = 1; j <= _r; j++) {
                    E_m += epsilons[j - 1] * (_E_t[j - 1] - _E_t[j]);
                }
                double log_E_m = std::log(E_m) + _epsilons_log_scale[m - 1];

                for (int i = 1; i <= m - 1; i++) {
                    double log_x_1 = (m - i - 1) * std::log(_E_t[_r]); // x_{m-i-1,r}
                    double log_var_x_m1 = std::log(_varX(i - 1, _r));
                    double log_x_2 = 2 * (i - 1) * std::log(_E_t[_r]);

                    std::vector<double> &epsilon = _epsilons_normalized[i - 1];
                    double tmp_sum = 0;
                    for (int j = 1; j <= _r; j++) {
                        double fac_1 = std::exp(log_x_2) * (_cov_t(_r, j - 1) - _cov_t(_r, j));
                        double fac_2 = std::exp(log_var_x_m1) * _E_t[_r] * (_E_t[j - 1] - _E_t[j]);
                        if (j == _r) {
                            fac_2 = std::exp(log_var_x_m1) * (_E_t[_r] * _E_t[j - 1] - _E_t2[j]);
                        }
                        tmp_sum += epsilon[j - 1] * (fac_1 + fac_2);
                    }
                    double log_e_sum = std::log(std::abs(tmp_sum)) + _epsilons_log_scale[i - 1];
                    int sgn_log_e_sum = tmp_sum < 0 ? -1 : 1;
                    log_parts[i - 1] = std::log(2) + log_E_m + log_x_1 + log_e_sum;
                    signs_parts[i - 1] = sgn_log_e_sum;
                }

                base::MathUtils::sumOfLog(log_parts, signs_parts, &log_part_1, &sign_1);
            }

            double log_var_x_m = std::log(_varX(m - 1, _r));
            double log_e_x_2 = (m - 1) * std::log(_E_t2[_r]);
            double tmp_sum_2 = 0;
            for (int j = 1; j <= _r; j++) {
                for (int k = 1; k <= _r; k++) {
                    double fac_1 = std::exp(log_e_x_2) * _cov_diff_t(j, k);
                    double fac_2 = std::exp(log_var_x_m) * (_E_t[j - 1] - _E_t[j]) * (_E_t[k - 1] - _E_t[k]);
                    tmp_sum_2 += epsilons[j - 1] * epsilons[k - 1] * (fac_1 + fac_2);
                }
            }
            double log_tmp_sum_2 = std::log(tmp_sum_2) + 2 * _epsilons_log_scale[m - 1];
            double log_part_2 = log_tmp_sum_2;
            int sign_2 = 1;


            base::MathUtils::sumOfLog({log_part_1, log_part_2}, {sign_1, sign_2}, &log_contributions[m - 1],
                                      &signs[m - 1]);
        }


        double log_zd_var = 0;
        int sign_zd = 1;
        base::MathUtils::sumOfLog(log_contributions, signs, &log_zd_var, &sign_zd);
        if (sign_zd < 0) {
            throw std::runtime_error("sign of variance is negative!!!");
        }
        return log_zd_var;

    }

    double PosteriorEstimator::_cov_t(int i, int j) {
        if (i == j) {
            double var_t = (_N - j + 1) * j / (double) ((_N + 2) * std::pow((double) _N + 1, 2));
            return var_t;
        } else {
            return 0;
        }

    }

    double PosteriorEstimator::_cov_diff_t(int i, int j) {
        double cov_diff = _cov_t(i - 1, j - 1) - _cov_t(i - 1, j) - _cov_t(i, j - 1) + _cov_t(i, j);
        return cov_diff;
    }


    void
    PosteriorEstimator::_computeCovTerm(DeadParticleSet &dead_points, LiveParticleSet &live_points, double *log_cov,
                                        int *sign_cov) {
        int m = dead_points.size() / _r;
        std::vector<double> log_parts(m, 0.0);
        std::vector<int> signs_parts(m, 1);

        double log_L_m = live_points.getLogAverageL();


        for (int i = 1; i <= m; i++) {
            double log_x_1 = (m - i) * std::log(_E_t[_r]);
            double log_var_x_m1 = std::log(_varX(i - 1, _r));
            double log_x_2 = 2 * (i - 1) * std::log(_E_t[_r]);

            std::vector<double> &epsilon = _epsilons_normalized[i - 1];
            double tmp_sum = 0;
            for (int j = 1; j <= _r; j++) {
                double fac_1 = std::exp(log_x_2) * (_cov_t(_r, j - 1) - _cov_t(_r, j));
                double fac_2 = std::exp(log_var_x_m1) * _E_t[_r] * (_E_t[j - 1] - _E_t[j]);
                if (j == _r) {
                    fac_2 = std::exp(log_var_x_m1) * (_E_t[_r] * _E_t[j - 1] - _E_t2[j]);
                }
                tmp_sum += epsilon[j - 1] * (fac_1 + fac_2);
            }
            double log_e_sum = std::log(std::abs(tmp_sum)) + _epsilons_log_scale[i - 1];
            int sgn_log_e_sum = tmp_sum < 0 ? -1 : 1;
            log_parts[i - 1] = std::log(2) + log_L_m + log_x_1 + log_e_sum;
            signs_parts[i - 1] = sgn_log_e_sum;
        }

        base::MathUtils::sumOfLog(log_parts, signs_parts, log_cov, sign_cov);
    }
}

