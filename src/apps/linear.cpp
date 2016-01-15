/**
 * git clone https://github.com/hou-dao/raman.git
 * ---
 * Writtened by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "deom.hpp"
#include "resp.hpp"

int main () {

    ifstream jsonFile("input.json");
    stringstream strStream;
    strStream << jsonFile.rdbuf();
    string jsonStr = strStream.str();
    string err;

    const Json json = Json::parse(jsonStr,err);
    if (!err.empty()) {
        printf ("Error in parsing input file: %s\n", err.c_str());
        return 0;
    }

    syst s(json["deom"]["syst"]);
    bath b(json["deom"]["bath"]);
    hidx h(json["deom"]["hidx"]);

    const double w_max = json["linear"]["wmax"].number_value();              
    const int    nt = json["linear"]["nt"].int_value();
    const double dt = json["linear"]["dt"].number_value();
    const double staticErr = json["linear"]["staticErr"].number_value();
    const int    nk = json["linear"]["nk"].int_value();
    const string dip0File = json["linear"]["dip0"].string_value();
    const string dip1File = json["linear"]["dip1"].string_value();
    cx_mat dipo0, dipo1;
    if (dipo0.load (dip0File, arma_ascii)) {
        dipo0.print("dipole 0");
    } else {
        printf("Fail to load dipo0!\n");
    }
    if (dipo1.load (dip1File, arma_ascii)) {
        dipo1.print("dipole 1");
    } else {
        printf("Fail to load dipo1!\n");
    }

    resp1st (w_max, nt, dt, staticErr, nk, dipo0, dipo1, s, b, h);

    return 0;
}
