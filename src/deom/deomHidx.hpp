/**
 * git clone https://github.com/hou-dao/raman.git
 * ---
 * Writtened by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#pragma once

#ifdef _MAP
#include <map>
#else
#include <unordered_map>
#endif

#include <cstdio>
#include <string>
#include "armadillo"
#include "json11.hpp"

using namespace std;
using namespace arma;
using namespace json11;

class hnod {

    public:

    int       rank;
    int       tier;
    cx_double gams;

    hnod (): rank(0), tier(0), gams(0) {};

    hnod (const hnod& rhs): rank(rhs.rank), tier(rhs.tier), gams(rhs.gams) {}

    hnod (const int _rank, const int _tier, const cx_double _gams): 
        rank(_rank), tier(_tier), gams (_gams) {};

    hnod& operator= (const hnod& rhs) {
        if (this != &rhs) {
            rank = rhs.rank;
            tier = rhs.tier;
            gams = rhs.gams;
        }
        return *this;
    };

   ~hnod () {rank = tier = -9527; gams = 0;};

};

typedef field<string> hkey;
#ifdef _MAP
typedef map<string,hnod> hmap;
#else
typedef unordered_map<string,hnod> hmap;
#endif

class hidx {

    public:

    int       nind;
    int       lmax;
    int       nmax;
    int       lddo;
    int       nddo;
    double    ferr;
    hkey      keys;
    hmap      maps;
    cx_rowvec expn;

    hidx (const Json& json) {

        const string expnFile = json["expnFile"].string_value();
        cx_mat expn_gam;
        if (expn_gam.load (expnFile, arma_ascii)) {
            expn = vectorise(expn_gam,1);
        } else {
            printf ("Error: expn_gam is not loaded!\n");
        }

        nind = expn.n_elem;
        lmax = json["lmax"].int_value();
        nmax = json["nmax"].int_value();
        ferr = json["ferr"].number_value();
        
        nddo = 1;
        lddo = 0;
        int ntot = static_cast<int>(get_nmax(nind,lmax));
        if (ntot < nmax) {
            nmax = ntot;
        }

        keys.set_size(nmax);
        keys(0) = string (nind,'0');
        for (int i=1; i<nmax; ++i) {
            keys(i) = string (nind,'!');
        }
        if (!maps.empty()) {
            printf ("Madan! maps is not empty!");
        }
        maps.emplace(keys(0),hnod());

        printf ("$InitHidx\n");
        printf ("nind = %d\n", nind);
        printf ("lddo = %d\n", lddo);
        printf ("nddo = %d\n", nddo);
        printf ("lmax = %d\n", lmax);
        printf ("nmax = %d\n", nmax);
        printf ("ferr = %g\n", ferr);
        printf ("$InitHidx\n\n");
    }

    hidx (const hidx& rhs): nind(rhs.nind), lmax(rhs.lmax),nmax(rhs.nmax),
        lddo (rhs.lddo), nddo(rhs.nddo), ferr(rhs.ferr), keys(rhs.keys), 
        maps(rhs.maps), expn(rhs.expn) {}

    hidx (const int _nind, const int _lmax, const int _nmax, const int _lddo, const int _nddo,
            const double _ferr, const hkey& _keys, const hmap& _maps, const cx_rowvec& _expn):
        nind(_nind), lmax(_lmax),nmax(_nmax),
        lddo (_lddo), nddo(_nddo), ferr(_ferr), keys(_keys), 
        maps(_maps), expn(_expn) {}

   ~hidx () {nind = nmax = lmax = nddo = lddo = -9527;}

    hidx& operator= (const hidx& rhs) {
        if (this != &rhs) {
            nind = rhs.nind;
            lmax = rhs.lmax;
            nmax = rhs.nmax;
            lddo = rhs.lddo;
            nddo = rhs.nddo;
            ferr = rhs.ferr;
            keys = rhs.keys;
            maps = rhs.maps;
            expn = rhs.expn;
        }
        return *this;
    }
    
    unsigned long get_nmax (const int K, const int L) const {
        unsigned long ntot = 1;
        for (int k=1; k<=K; ++k) {
            ntot *= L+k;
            ntot /= k;
            if (ntot > 300000) {
                printf ("Be careful! too many elements");
            }
        }
        return ntot;
    }
};
