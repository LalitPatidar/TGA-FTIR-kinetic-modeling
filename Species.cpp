#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include "Species.h"
#include <map>

Species::Species(std::string &line1, std::string &line2, std::string &line3, std::string &line4)
{
    std::stringstream ss(line1);
    std::vector <std::string> vline1;
    std::string item;

    while(getline(ss,item,' '))
    {
        if (!item.empty())
        {
            vline1.push_back(item);
        }
    }
    
    name = vline1[0];
    int nC = atoi(vline1[2].substr(0,vline1[2].size()-1).c_str());
    int nH = atoi(vline1[3].substr(0,vline1[3].size()-1).c_str());
    int nN = atoi(vline1[4].substr(0,vline1[4].size()-1).c_str());
    int nO = atoi(vline1[5].substr(0,vline1[5].size()-1).c_str());

    composition.insert(std::make_pair("C",nC));
    composition.insert(std::make_pair("H",nH));
    composition.insert(std::make_pair("N",nN));
    composition.insert(std::make_pair("O",nO));
    
    molecular_weight = 12.0107*nC + 1.0079*nH + 14.0067*nN + 15.9994*nO;
}


