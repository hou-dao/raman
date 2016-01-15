/**
 * git clone https://github.com/hou-dao/raman.git
 * ---
 * Writtened by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#pragma once

#include "deom.hpp"

void resp1st (const double w_max, const int nt, const double dt,
              const double staticErr, const int nk,
              const cx_mat& dipo0, const cx_mat& dipo1, 
              const syst& s, const bath& b, const hidx& h);

void resp2nd (const double w1_max, const double w2_max, 
              const int nt1, const int nt2, const double dt,
              const double staticErr, const int nk,
              const cx_mat& dipo0, const cx_mat& dipo1, cx_mat& dipo2, 
              const char sch_hei, const syst& s, const bath& b, const hidx& h);
