//
// Created by jan on 12/10/18.
//

#ifndef LFNS_PARSERREADER_H
#define LFNS_PARSERREADER_H

#include <map>
#include "TxtFileReader.h"

namespace io {
    const std::string PARAMETER_CAPTION = "Parameters:";
    const std::string SPECIES_CAPTION = "Species:";
    const std::string RANDOM_NUMBER_CAPTION = "Random numbers:";

    enum class RNDNBR_TYPE {
        NORMAL, UNIFORM, UNIFORM_INT, POISSON
    };
    //TODO think of more elegant way...
    static std::map<RNDNBR_TYPE, std::string> RNDNBR_NAME = {{ RNDNBR_TYPE::NORMAL, "Normal"},
                                                                {RNDNBR_TYPE::UNIFORM, "Uniform"},
                                                                 {RNDNBR_TYPE::UNIFORM_INT, "Uniform_Int"},
                                                                 {RNDNBR_TYPE::POISSON, "Poisson"}};

    class ParserReader : public TxtFileReader {
    public:
        ParserReader(std::string parser_file_name);

        virtual ~ParserReader();

        bool parametersDefined() const;

        bool speciesDefined() const;

        bool randomNumbersDefined() const;

        std::vector<std::string> readParameterNames();

        std::vector<std::string> readSpeciesNames();

        std::vector<std::string> readRdn(RNDNBR_TYPE type);
        std::vector<std::string> readNormalRdn();

        std::vector<std::string> readUniformRdn();

        std::vector<std::string> readUniformIntRdn();
        std::vector<std::string> readPoissonRdn();

        std::vector<std::pair<double, double> > readNormalParams();

        std::vector<std::pair<double, double> > readUniformParams();

        std::vector<std::pair<int, int> > readUniformIntParams();

        std::vector<double> readPoissonParams();

    };
}


#endif //LFNS_PARSERREADER_H
