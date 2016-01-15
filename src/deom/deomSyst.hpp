/**
 * git clone https://github.com/hou-dao/raman.git
 * ---
 * Writtened by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#pragma once

#include <string>
#include "armadillo"
#include "json11.hpp"

using namespace std;
using namespace arma;
using namespace json11;

class syst {
    
    public:
    
    cx_mat  ham1;
    cx_cube qmd1;
    
    syst (const Json& json) {
        const string hamsFile = json["hamsFile"].string_value();
        const string qmdsFile = json["qmdsFile"].string_value();
        printf ("$InitSyst\n");
        if (ham1.load (hamsFile, arma_ascii)) {
            ham1.print("ham1");
        } else {
            printf ("Fail to load ham1!\n");
        }
        if (qmd1.load (qmdsFile, arma_ascii)) {
            qmd1.print("qmd1");
        } else {
            printf ("Fail to load qmd1!\n");
        }
        printf ("$InitSyst\n\n");
    }
    
    syst (const cx_mat& _h, const cx_cube& _q): ham1(_h), qmd1(_q) {}
        
   ~syst () {}
};
