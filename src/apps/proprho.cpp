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

    deom d(json["deom"]);

    const int inistate = json["proprho"]["inistate"].int_value();
    const int    nt = json["proprho"]["nt"].int_value();
    const double dt = json["proprho"]["dt"].number_value();
    const int    nk = json["proprho"]["nk"].int_value();

    cx_cube ddos = zeros<cx_cube>(size(d.ddos1));
    ddos(inistate,inistate,0) = 1;
    d.propagation (ddos, nt, dt, nk);

    return 0;
}
