#ifndef SPECIES_H
#define SPECIES_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <map>

class Species {
    public:
        Species(std::string &line1, std::string &line2, std::string &line3, std::string &line4);
        std::string name;
        std::map<std::string, int> composition;
        double molecular_weight;
        double radius;
    private:
};

#endif
