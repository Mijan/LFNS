//
// Created by jan on 31/10/18.
//

#ifndef LFNS_PARTICLEFILTERSETTINGS_H
#define LFNS_PARTICLEFILTERSETTINGS_H


#include <string>
#include <map>
#include <vector>
#include <iostream>

namespace particle_filter {
    class ParticleFilterSettings {
    public:
        int num_used_trajectories = 1e5;
        std::map<std::string, std::string> data_files;
        std::map<std::string, std::string> time_files;

        bool use_premature_cancelation = false;
        int H = 200;

        std::string provided_parameter_file = "";
        std::vector<std::string> provided_parameter_names = {};

        void print(std::ostream &stream) {


            stream << "H:\t" << H << std::endl;
            stream << "Data used:" << std::endl;
            std::map<std::string, std::string>::iterator it;
            for (it = data_files.begin(); it != data_files.end(); it++) {
                stream << "Experiment:\t" << it->first << "\t\tData file: " << it->second
                       << "\ttimes file: " << time_files[it->first] << std::endl;
            }
            stream << "Nnumber of trajectories used for each experiment: " << num_used_trajectories << std::endl;
            stream << "\n\n" << std::endl;


            if (provided_parameter_file.size() > 0) {
                stream << "\nA file with already presampled parameters was provided the parameters" << std::endl;
                for (std::string &param_name : provided_parameter_names) stream << param_name << " ";
                stream << "\nwill be read from the file" << std::endl;
                stream << provided_parameter_file << std::endl;
            }
        }

        std::vector<std::string> param_names;
    };
}


#endif //LFNS_PARTICLEFILTERSETTINGS_H
