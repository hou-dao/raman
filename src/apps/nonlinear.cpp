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

    const double w1_max = json["nonlinear"]["w1max"].number_value();
    const double w2_max = json["nonlinear"]["w2max"].number_value();
    const int    nt1 = json["nonlinear"]["nt1"].int_value();
    const int    nt2 = json["nonlinear"]["nt1"].int_value();
    const double dt = json["nonlinear"]["dt"].number_value();
    const double staticErr = json["nonlinear"]["staticErr"].number_value();
    const int    nk = json["nonlinear"]["nk"].int_value();
    const string dip0File = json["nonlinear"]["dip0"].string_value();
    const string dip1File = json["nonlinear"]["dip1"].string_value();
    const string dip2File = json["nonlinear"]["dip1"].string_value();
    const string sch_hei = json["nonlinear"]["sch_hei"].string_value();
    
    cx_mat dipo0, dipo1, dipo2;
    if (!dipo0.load (dip0File, arma_ascii)) {
        printf ("Fail to load dipo0!\n");
    } else {
        dipo0.print ("dipole 0");    
    }
    if (!dipo1.load (dip1File, arma_ascii)) {
        printf ("Fail to load dipo1!\n");
    } else {
        dipo1.print ("dipole 1");    
    }
    if (!dipo2.load (dip2File, arma_ascii)) {
        printf ("Fail to load dipo2!\n");
    } else {
        dipo2.print ("dipole 2");
    }

    resp2nd (w1_max, w2_max, nt1, nt2, dt, staticErr, nk,
             dipo0, dipo1, dipo2, sch_hei[0], s, b, h);

    return 0;
}
