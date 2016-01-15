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
#include "deomSyst.hpp"
#include "deomBath.hpp"
#include "deomHidx.hpp"

using namespace std;
using namespace arma;
using namespace json11;

class deom: public syst, public bath, public hidx {

    public:

    // Syst and Bath
    int nsys;
    int nmod;
    int nper;
    
    // Action 
    cx_cube qddo;
    cx_cube ddoq;

    // Auxiliary for rk4
    cx_cube ddos1;
    cx_cube ddos2;
    cx_cube ddos3;

    deom (const Json& json): syst (json["syst"]), bath (json["bath"]), hidx (json["hidx"]) {
        nsys = ham1.n_rows;
        nmod = qmd1.n_slices;
        nper = nind/nmod;
        qddo.copy_size(qmd1);
        ddoq.copy_size(qmd1);
        ddos1.set_size(nsys,nsys,nmax);
        ddos2.set_size(nsys,nsys,nmax);
        ddos3.set_size(nsys,nsys,nmax);
    }

    deom (const syst& s, const bath& b, const hidx& h): syst (s), bath (b), hidx (h) {
        nsys = ham1.n_rows;
        nmod = qmd1.n_slices;
        nper = nind/nmod;
        qddo.copy_size(qmd1);
        ddoq.copy_size(qmd1);
        ddos1.set_size(nsys,nsys,nmax);
        ddos2.set_size(nsys,nsys,nmax);
        ddos3.set_size(nsys,nsys,nmax);
    }
    
    deom (const deom& rhs): syst(rhs.ham1,rhs.qmd1), 
          bath(rhs.temperature, rhs.coef_lft, rhs.coef_rht, rhs.coef_abs, rhs.expn_gam, rhs.delt_res), 
          hidx(rhs.nind, rhs.lmax, rhs.nmax, rhs.lddo, rhs.nddo, 
               rhs.ferr, rhs.keys, rhs.maps, rhs.expn) {
        nsys = ham1.n_rows;
        nmod = qmd1.n_slices;
        nper = nind/nmod;
        qddo.copy_size(qmd1);
        ddoq.copy_size(qmd1);
        ddos1.set_size(nsys,nsys,nmax);
        ddos2.set_size(nsys,nsys,nmax);
        ddos3.set_size(nsys,nsys,nmax);              
    }

   ~deom () {}

    void rem (cx_cube& d_ddos, const cx_cube& ddos, const double t) {
        rem_sch (d_ddos, ddos, t);
    }

    void rem (cx_cube& d_ddos, const cx_cube& ddos, const double t, const char sch_hei) {
        if (sch_hei == 's') {
            rem_sch (d_ddos, ddos, t);
        } else if (sch_hei == 'h') {
            rem_hei (d_ddos, ddos, t);
        } else {
            printf("sch_hei is invalid!\n");
        }
    }

    void rem_sch (cx_cube& d_ddos, const cx_cube& ddos, const double t);

    void rem_hei (cx_cube& d_ddos, const cx_cube& ddos, const double t);

    void propagation (cx_cube& ddos, const double dt=0.005, const int nt=1000, const int nk=10);

    void equilibrium (cx_cube& ddos, const double dt=0.005, const double err=2.e-8, const int nk=10);

    inline bool is_valid (const cx_mat& ddo) const {
        return any(abs(vectorise(ddo))>ferr);
    }

    void filter (cx_cube& ddos) {
        int n = 1;
        int l = 0;
        for (int iddo=1; iddo<nddo; ++iddo) {
            auto iter = maps.find(keys(iddo));
		    if (is_valid(ddos.slice(iddo))) {
			    if (n != iddo) {
                    iter->second.rank = n;
                    keys(n) = keys(iddo);
                    ddos.slice(n) = ddos.slice(iddo);
			    }
			    ++n;
		        l = l>(iter->second.tier)?l:(iter->second.tier);
		    } else {
			    maps.erase(iter);
            }
	    } 
        lddo = l;
        nddo = n;
    }

    template<typename... Tc>
    void rk4 (cx_cube& ddos, const double t, const double dt, const Tc&... args) {
    
        const double dt2 = dt*0.5;
        const double dt6 = dt/6.0;
    
        cx_cube ddos1(size(ddos));
        cx_cube ddos2(size(ddos));
        cx_cube ddos3(size(ddos));
    
        // K1
        const int nddo0 = nddo;
        rem (ddos1,ddos,t,args...);
        ddos3.slices(0,nddo0-1) = ddos.slices(0,nddo0-1)+ddos1.slices(0,nddo0-1)*dt2;
        if (nddo > nddo0) {
            ddos3.slices(nddo0,nddo-1) = ddos1.slices(nddo0,nddo-1)*dt2;
        }
        // K2
        const int nddo1 = nddo;
        rem (ddos2,ddos3,t+0.5*dt,args...);
        ddos1.slices(0,nddo1-1) += ddos2.slices(0,nddo1-1)*2.0;
        if (nddo > nddo1) {
            ddos1.slices(nddo1,nddo-1) = ddos2.slices(nddo1,nddo-1)*2.0;
        }
        ddos3.slices(0,nddo0-1) = ddos.slices(0,nddo0-1)+ddos2.slices(0,nddo0-1)*dt2;
        if (nddo > nddo0) {
            ddos3.slices(nddo0,nddo-1) = ddos2.slices(nddo0,nddo-1)*dt2;
        }
        // K3
        const int nddo2 = nddo;
        rem (ddos2,ddos3,t+0.5*dt,args...);
        ddos1.slices(0,nddo2-1) += ddos2.slices(0,nddo2-1)*2.0;
        if (nddo > nddo2) {
            ddos1.slices(nddo2,nddo-1) = ddos2.slices(nddo2,nddo-1)*2.0;
        }
        ddos3.slices(0,nddo0-1) = ddos.slices(0,nddo0-1)+ddos2.slices(0,nddo0-1)*dt;
        if (nddo > nddo0) {
            ddos3.slices(nddo0,nddo-1) = ddos2.slices(nddo0,nddo-1)*dt;
        }
        // K4
        const int nddo3 = nddo;
        rem (ddos2,ddos3,t+dt,args...);
        ddos1.slices(0,nddo3-1) += ddos2.slices(0,nddo3-1);
        if (nddo > nddo3) {
            ddos1.slices(nddo3,nddo-1) = ddos2.slices(nddo3,nddo-1);
        }
        ddos.slices(0,nddo0-1) += ddos1.slices(0,nddo0-1)*dt6;
        if (nddo > nddo0) {
            ddos.slices(nddo0,nddo-1) = ddos1.slices(nddo0,nddo-1)*dt6;
        }
        filter (ddos);
    }
};
