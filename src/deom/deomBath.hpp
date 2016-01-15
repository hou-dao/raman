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
#include "deomConst.hpp"

using namespace std;
using namespace arma;
using namespace json11;

class bath {

    public:

    double temperature;
    cx_mat coef_lft;
    cx_mat coef_rht;
    mat    coef_abs;
    cx_mat expn_gam;
    vec    delt_res;

    bath (const Json& json) {

        const string etalFile = json["etalFile"].string_value();
        const string etarFile = json["etarFile"].string_value();
        const string etaaFile = json["etaaFile"].string_value();
        const string expnFile = json["expnFile"].string_value();
        const string delrFile = json["delrFile"].string_value();

        temperature = json["temp"].number_value();

        printf ("$InitBath\n");
        if (coef_lft.load (etalFile, arma_ascii)) {
            coef_lft.print("coef_lft");
        } else {
            printf ("coef_lft is not loaded!\n");
        }
        if (coef_rht.load (etarFile, arma_ascii)) {
            coef_rht.print("coef_rht");
        } else {
            printf ("coef_rht is not loaded!\n");
        }
        if (coef_abs.load (etaaFile, arma_ascii)) {
            coef_abs.print("coef_abs");
        } else {
            printf ("coef_abs is not loaded!\n");
        }
        if (expn_gam.load (expnFile, arma_ascii)) {
            expn_gam.print("expn_gam");
        } else {
            printf ("expn_gam is not loaded!\n");
        }
        if (delt_res.load (delrFile, arma_ascii)) {
            delt_res.print("delt_res");
        } else {
            printf ("delt_res is not loaded!\n");
        }
        printf ("$InitBath\n\n");
    }

    bath (const bath& rhs): temperature(rhs.temperature),
          coef_lft(rhs.coef_lft), coef_rht(rhs.coef_rht), 
          coef_abs(rhs.coef_abs), expn_gam(rhs.expn_gam), 
          delt_res(rhs.delt_res) {};

    bath (const double temp, const cx_mat& etal, const cx_mat& etar, const mat& etaa, 
          const cx_mat& expn, const vec& delr): temperature(temp),
          coef_lft(etal), coef_rht(etar), coef_abs(etaa),
          expn_gam(expn), delt_res(delr) {};

    bath& operator= (const bath& rhs) {
        if (this != &rhs) {
            temperature = rhs.temperature;
            coef_lft = rhs.coef_lft;
            coef_rht = rhs.coef_rht;
            coef_abs = rhs.coef_abs;
            expn_gam = rhs.expn_gam;
            delt_res = rhs.delt_res;
        }
        return *this;
    }

   ~bath () {};

};
