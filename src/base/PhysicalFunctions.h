//
// Created by jan on 5/16/23.
//

#ifndef LFNS_PHYSICALFUNCTIONS_H
#define LFNS_PHYSICALFUNCTIONS_H


namespace base {
    class PhysicalFunctions {
    public:
        static double averageCylindricLightIntensity(double X, double alpha, double R, double beta, double I0);

        static double
        cylindricLightIntensityIntegrand(double phi, double z, double X, double alpha, double R, double beta);

        static double cylindricIntegratePhi(double z, double X, double alpha, double R, double beta);

    };
}

#endif //LFNS_PHYSICALFUNCTIONS_H
