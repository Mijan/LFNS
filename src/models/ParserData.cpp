//
// Created by jan on 12/10/18.
//

#include <stdexcept>
#include <sstream>
#include "ParserData.h"
#include "../io/ParserReader.h"

namespace models {
    ParserData::ParserData(std::vector<std::string> parameter_names_, std::vector<std::string> species_names_)
            : BaseData(parameter_names_, species_names_), _normal_random_nbr(), _uniform_random_nbr(),
              _uniform_int_random_nbr(), _poisson_random_nbr(), _normal_random_params(), _uniform_random_params(),
              _uniform_int_random_params(), _poisson_random_params() {}

    ParserData::ParserData(std::string file_name) : BaseData(), _normal_random_nbr(), _uniform_random_nbr(),
                                                    _uniform_int_random_nbr(), _poisson_random_nbr(),
                                                    _normal_random_params(), _uniform_random_params(),
                                                    _uniform_int_random_params(), _poisson_random_params() {
        try {
            io::ParserReader reader(file_name);
            if (reader.parametersDefined()) { parameter_names = reader.readParameterNames(); }
            if (reader.speciesDefined()) { species_names = reader.readSpeciesNames(); }
            if (reader.randomNumbersDefined()) {
                _normal_random_nbr = reader.readNormalRdn();
                if (!_normal_random_nbr.empty()) { _normal_random_params = reader.readNormalParams(); }
                _uniform_random_nbr = reader.readUniformRdn();
                if (!_uniform_random_nbr.empty()) { _uniform_random_params = reader.readUniformParams(); }
                _uniform_int_random_nbr = reader.readUniformIntRdn();
                if (!_uniform_int_random_nbr.empty()) { _uniform_int_random_params = reader.readUniformIntParams(); }
                _poisson_random_nbr = reader.readPoissonRdn();
                if (!_poisson_random_nbr.empty()) { _poisson_random_params = reader.readPoissonParams(); }
                if (_normal_random_nbr.size() != _normal_random_params.size()) {
                    std::stringstream ss;
                    ss << "Failed to define random numbers. Number of normal random numbers defined: "
                       << _normal_random_nbr.size() << ", but " << _normal_random_params.size()
                       << " parameters defined!" << std::endl;
                    throw std::runtime_error(ss.str());
                }
                if (_uniform_random_nbr.size() != _uniform_random_params.size()) {
                    std::stringstream ss;
                    ss << "Failed to define random numbers. Number of uniform random numbers defined: "
                       << _uniform_random_nbr.size() << ", but " << _uniform_random_params.size()
                       << " parameters defined!" << std::endl;
                    throw std::runtime_error(ss.str());
                }
                if (_uniform_int_random_nbr.size() != _uniform_int_random_params.size()) {
                    std::stringstream ss;
                    ss << "Failed to define random numbers. Number of uniform int random numbers defined: "
                       << _uniform_int_random_nbr.size() << ", but " << _uniform_int_random_params.size()
                       << " parameters defined!" << std::endl;
                    throw std::runtime_error(ss.str());
                }
                if (_poisson_random_nbr.size() != _poisson_random_params.size()) {
                    std::stringstream ss;
                    ss << "Failed to define random numbers. Number of uniform int random numbers defined: "
                       << _poisson_random_nbr.size() << ", but " << _poisson_random_params.size()
                       << " parameters defined!" << std::endl;
                    throw std::runtime_error(ss.str());
                }
            }
        } catch (const std::exception &e) {
            std::stringstream os;
            os << "Failed to crate Parser Data: " << "\n\t" << e.what() << std::endl;
            throw std::runtime_error(os.str());
        }
    }


    ParserData::ParserData(const ParserData &data) : BaseData(data),
                                                     _normal_random_nbr(data._normal_random_nbr),
                                                     _uniform_random_nbr(data._uniform_random_nbr),
                                                     _uniform_int_random_nbr(data._uniform_int_random_nbr),
                                                     _poisson_random_nbr(data._poisson_random_nbr),
                                                     _normal_random_params(data._normal_random_params),
                                                     _uniform_random_params(data._uniform_random_params),
                                                     _uniform_int_random_params(data._uniform_int_random_params),
                                                     _poisson_random_params(data._poisson_random_params) {}

    ParserData::~ParserData() {}

    const std::vector<std::string> &ParserData::getNormalRandomNbrsName() const {
        return _normal_random_nbr;
    }

    const std::vector<std::string> &ParserData::getUniformRandomNbrsName() const {
        return _uniform_random_nbr;
    }

    const std::vector<std::string> &ParserData::getUniformIntRandomNbrsName() const {
        return _uniform_int_random_nbr;
    }

    const std::vector<std::string> &ParserData::getPoissonRandomNbrsName() const {
        return _poisson_random_nbr;
    }

    const std::vector<std::pair<double, double> > &ParserData::getNormalRandomParams() const {
        return _normal_random_params;
    }

    const std::vector<std::pair<double, double> > &ParserData::getUniformRandomParams() const {
        return _uniform_random_params;
    }

    const std::vector<std::pair<int, int> > &ParserData::getUniformIntRandomParams() const {
        return _uniform_int_random_params;
    }

    const std::vector<double> &ParserData::getPoissonRandomParams() const {
        return _poisson_random_params;
    }

    std::size_t ParserData::getNumNormalNumbers() { return _normal_random_nbr.size(); }

    std::size_t ParserData::getNumUniformNumbers() { return _uniform_random_nbr.size(); }

    std::size_t ParserData::getNumUniformIntNumbers() { return _uniform_int_random_nbr.size(); }

    std::size_t ParserData::getNumPoissonNumbers() { return _poisson_random_nbr.size(); }

}