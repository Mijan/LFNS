//
// Created by jan on 24/10/18.
//

#ifndef LFNS_CONFIGFILEINTERPRETER_H
#define LFNS_CONFIGFILEINTERPRETER_H

#include <string>
#include <vector>
#include <map>
#include "XmlFileReader.h"
#include "XmlPropertyMap.h"
#include "../base/Utils.h"

namespace io {

    class ConfigFileException : public std::exception {
    public:
        ConfigFileException(std::string e) : _error_message(e) {};

        const char *what() const throw() { return _error_message.c_str(); }

    private:
        std::string _error_message;
    };

    class ConfigFileInterpreter {
    public:
        ConfigFileInterpreter(std::string instance_file_name);

        virtual ~ConfigFileInterpreter();

        std::string getModelFileName();

        std::string getMeasurementModelFile();

        std::string getInitialConditionsFile();

        std::string getModelType();

        std::map<std::string, std::string> getParameterScales();

        std::map<std::string, std::pair<double, double> > getParameterBounds();

        std::map<std::string, double> getFixedParameters();

        std::map<std::string, std::string> getDataFiles(std::vector<std::string> experiments);

        std::map<std::string, std::string> getTimesFiles(std::vector<std::string> experiments);

        std::vector<std::string> getExperimentsForLFNS();

        std::vector<std::string> getExperimentsForSimulations();

        std::vector<std::string> getExperimentsForEvaluateLikelihood();

        int getNForLFNS();

        int getRForLFNS();

        int getHForLFNS();

        int getHForEvaluateLikelihood();

        int getDPGMMItForLFNS();

        double getEpsilonForLFNS();

        std::vector<double> getPulsePeriods(std::string experiment_name);

        std::vector<double> getPulseStrengths(std::string experiment_name);

        std::vector<double> getPulseDurations(std::string experiment_name);

        std::vector<int> getNumPulse(std::string experiment_name);

        std::vector<std::string> getPulseInputNames(std::string experiment_name);

        std::vector<double> getStartingTimes(std::string experiment_name);

        std::vector<double> getParamForSimulation();

        std::vector<double> getParamForEvaluateLikelihood();

        int getNForSimulation();

        std::string getParameterFileforSimulation();

        std::string getParameterFileforEvaluateLikelihood();

        double getInitialTimeForSimulation();

        double getFinalTimeForSimulation();

        double getIntervalForSimulation();

        bool detSpeciesProvided();

        bool stochSpeciesProvided();

        std::vector<std::string> getDetSpecies();

        std::vector<std::string> getStochSpecies();

        double getMinStepSize();

        double getRelTol();

        double getAbsTol();

        int getMaxErrorFails();

        int getMaxNumberSteps();

        bool parametersProvided();

        std::string gerProvidedParametersFile();

        std::vector<std::string> getProvidedParameters();

        std::vector<std::string> getStepInputNames(std::string basicString);

        std::vector<std::vector<double>> getStepInputTimes(std::string experiment);

        std::vector<std::vector<double>> getStepInputStrengths(std::string experiment);

    protected:
        XmlFileReader _reader;

        std::vector<std::string> _getPulseInputValues(std::string experiment_name, std::string field_name);
        std::vector<std::string> _getStepInputValues(std::string experiment_name, std::string field_name);
        std::vector<std::string> _getInputValues(std::string input_key, std::string experiment_name, std::string field_name, bool split_values=true);


    };
}


#endif //LFNS_CONFIGFILEINTERPRETER_H
