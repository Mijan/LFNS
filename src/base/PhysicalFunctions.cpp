//
// Created by jan on 5/16/23.
//

#include <cmath>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "PhysicalFunctions.h"

namespace base {
    double PhysicalFunctions::cylindricLightIntensityIntegrand(double phi, double z, double X, double alpha, double R,
                                                               double beta) {
        double tmp = (R - z) * cos(phi) + sqrt(R * R - (R - z) * sin(phi) * sin(phi)) + beta;
        return exp(-alpha * X * tmp);
    }


    double PhysicalFunctions::cylindricIntegratePhi(double z, double X, double alpha, double R, double beta) {
        auto integrandFunc = [&](double phi) {
            return cylindricLightIntensityIntegrand(phi, z, X, alpha, R, beta);
        };

        boost::math::quadrature::gauss_kronrod<double, 15> integrator;
        double result = integrator.integrate(integrandFunc, 0, (double) M_PI);

        return result;
    }


    double PhysicalFunctions::averageCylindricLightIntensity(double X, double alpha, double R, double beta, double I0) {
        auto integrandFunc = [&](double z) {
            return cylindricIntegratePhi(z, X, alpha, R, beta);
        };

        boost::math::quadrature::gauss_kronrod<double, 15> integrator;
        double result = integrator.integrate(integrandFunc, 0, R);

        return I0 * result / ((double) M_PI * R);
    }
}