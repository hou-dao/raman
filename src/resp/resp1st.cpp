/**
 * git clone https://github.com/hou-dao/raman.git
 * ---
 * Writtened by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cmath>
#include <sstream>
#include "resp.hpp"

void resp1st (const double w_max, const int nt, const double dt,
              const double staticErr, const int nk,
              const cx_mat& dipo0, const cx_mat& dipo1, 
              const syst& s, const bath& b, const hidx& h) {

    const double dw1 = w_max/nt;
    const double dt1 = 2.0*deom_pi/w_max;
    const int    mt  = floor(dt1/dt);
    const double dt1_res = dt1-dt*mt; 

    deom d(s,b,h);

    const mat& exph= expmat(-real(d.ham1)/d.temperature);
    cx_cube rhot = zeros<cx_cube>(size(d.ddos1));
    rhot.slice(0).set_real(exph/trace(exph));
    d.equilibrium (rhot,dt,staticErr,nk);
    for (int iado=0; iado<d.nddo; ++iado) {
        rhot.slice(iado) = dipo0*rhot.slice(iado)-rhot.slice(iado)*dipo0;
    }

    cx_vec ft = zeros<cx_vec>(nt); 
    
    for (int it=0; it<nt; ++it) {
        ft(it) = trace(dipo1*rhot.slice(0));
        printf ("In sch-pic: it=%d, nddo=%d, lddo=%d\n", it, d.nddo, d.lddo);
        double t1 = it*dt1;
        for (int jt=0; jt<mt; ++jt) {
            d.rk4 (rhot,t1,dt);
            t1 += dt;
        }
        d.rk4 (rhot,t1,dt1_res);
    }

    // write time-domain signal
    const vec& ft_t1 = linspace(0.0,dt1*(nt-1)/deom_fs2unit,nt);
    const vec& ft_re = real(ft);
    const vec& ft_im = imag(ft);
    ft_t1.save("resp1st.t1",raw_ascii);
    ft_re.save("resp1st_re.t",raw_ascii);
    ft_im.save("resp1st_im.t",raw_ascii);
    
    // 1D FFT
    ft(0) *= 0.5;
    const cx_vec& fw = ifft(ft)*nt*dt1;
    
    // write freq-domain signal
    const vec& fw_w1 = linspace(0.0,dw1*(nt/2-1)/deom_cm2unit,nt/2);
    const vec& fw_re = real(fw.rows(0,nt/2-1));
    const vec& fw_im = imag(fw.rows(0,nt/2-1));
    fw_w1.save("resp1st.w1",raw_ascii);
    fw_re.save("resp1st_re.w",raw_ascii);
    fw_im.save("resp1st_im.w",raw_ascii);
}
