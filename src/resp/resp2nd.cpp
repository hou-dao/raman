/**
 * git clone https://github.com/hou-dao/raman.git
 * ---
 * Writtened by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cmath>
#include <sstream>
#include "resp.hpp"

void resp2nd (const double w1_max, const double w2_max,
              const int nt1, const int nt2, const double dt,
              const double staticErr, const int nk,
              const cx_mat& dipo0, const cx_mat& dipo1, cx_mat& dipo2, 
              const char sch_hei, const syst& s, const bath &b, const hidx& h) {

    const double dw1 = w1_max/nt1;
    const double dw2 = w2_max/nt2;
    const double dt1 = 2.0*deom_pi/w1_max;
    const double dt2 = 2.0*deom_pi/w2_max;
    const int    mt1 = floor(dt1/dt);
    const int    mt2 = floor(dt2/dt);
    const double dt1_res = dt1-dt*mt1;    
    const double dt2_res = dt2-dt*mt2;    
    
    deom d1(s,b,h);

    const mat& exph = expmat(-real(d1.ham1)/d1.temperature);
    cx_cube rho_t1 = zeros<cx_cube>(size(d1.ddos1));
    rho_t1.slice(0).set_real(exph/trace(exph));
    d1.equilibrium (rho_t1,dt,staticErr,nk);
    for (int iado=0; iado<d1.nddo; ++iado) {
        rho_t1.slice(iado) = dipo0*rho_t1.slice(iado)-rho_t1.slice(iado)*dipo0;
    }

    cx_mat ft = zeros<cx_mat>(nt2,nt1);

    if (sch_hei == 's') { // sch-pic
        
        for (int it1=0; it1<nt1; ++it1) {
            
            deom d2(d1);
            cx_cube rho_t2 = zeros<cx_cube>(size(d2.ddos1));
            for (int iado=0; iado<d1.nddo; ++iado) {
                rho_t2.slice(iado) = dipo1*rho_t1.slice(iado)-rho_t1.slice(iado)*dipo1;
            }

            for (int it2=0; it2<nt2; ++it2) {
                double t2 = it2*dt2;                    
                ft(it2,it1) = trace(dipo2*rho_t2.slice(0));
                printf("In sch-pic: it2=%d, nddo=%d, lddo=%d\n", it2, d2.nddo, d2.lddo);
                for (int jt2=0; jt2<mt2; ++jt2) {
                    d2.rk4 (rho_t2,t2,dt);
                    t2 += dt;
                }
                d2.rk4 (rho_t2,t2,dt2_res);
                t2 += dt2_res;
            }
            
            printf("In sch-pic: it1=%d, nddo=%d, lddo=%d\n", it1, d1.nddo, d1.lddo);
            double t1 = it1*dt1;
            for (int jt1=0; jt1<mt1; ++jt1) {
                d1.rk4 (rho_t1,t1,dt);
                t1 += dt;
            }
            d1.rk4 (rho_t1,t1,dt1_res);
            t1 += dt1_res;
        }
        
    } else if (sch_hei == 'h') { // hei-pic
    
        deom d2(s,b,h);
        cx_cube opr_t2 = zeros<cx_cube>(size(d2.ddos1));
        opr_t2.slice(0) = dipo2;
        for (int it2=0; it2<nt2; ++it2) {
            
            const hkey& keys = d2.keys.rows(0,d2.nddo-1);
            const cx_cube& oprs = opr_t2.slices(0,d2.nddo-1);
            stringstream ss1, ss2;
            ss1 << "key_t2_" << it2 << ".tmp";
            ss2 << "opr_t2_" << it2 << ".tmp";
            keys.save(ss1.str(),arma_binary);
            oprs.save(ss2.str(),arma_binary);

            printf("In hei-pic: it2=%d, nddo=%d, lddo=%d\n", it2, d2.nddo, d2.lddo);

            double t2 = it2*dt2;
            for (int jt2=0; jt2<mt2; ++jt2) {
                d2.rk4 (opr_t2,t2,dt,'h');
                t2 += dt;
            }
            d2.rk4 (opr_t2,t2,dt2_res,'h');
        }
        
        for (int it1=0; it1<nt1; ++it1) {

            cx_cube rho_t2 = zeros<cx_cube>(d1.nsys,d1.nsys,d1.nddo);
            for (int iado=0; iado<d1.nddo; ++iado) {
                rho_t2.slice(iado) = dipo1*rho_t1.slice(iado)-rho_t1.slice(iado)*dipo1;
            }
            
            for (int it2=0; it2<nt2; ++it2) {
                
                stringstream ss1, ss2;
                hkey    keys;
                cx_cube oprs;
                ss1 << "key_t2_" << it2 << ".tmp";
                ss2 << "opr_t2_" << it2 << ".tmp";
                keys.load(ss1.str(),arma_binary);
                oprs.load(ss2.str(),arma_binary);

                const int nddo_t2 = keys.n_rows;
                cx_double ctmp = 0.0;   
                for (int iado=0; iado<nddo_t2; ++iado) {
                    auto iter = d1.maps.find(keys(iado));
                    if (iter != d1.maps.end()) {
                        const int jado = iter->second.rank;
                        ctmp += trace(oprs.slice(iado)*rho_t2.slice(jado));
                    }
                }
                ft(it2,it1) = ctmp; 
            }
            
            printf("In sch-pic: it1=%d, nddo=%d, lddo=%d\n", it1, d1.nddo, d1.lddo);

            double t1 = it1*dt1;
            for (int jt1=0; jt1<mt1; ++jt1) {
                d1.rk4 (rho_t1,t1,dt);
                t1 += dt;
            }
            d1.rk4 (rho_t1,t1,dt1_res);
        }
       
        // clean tmp files
        for (int it2=0; it2<nt2; ++it2) {
            stringstream ss1, ss2;
            ss1 << "key_t2_" << it2 << ".tmp";
            ss2 << "opr_t2_" << it2 << ".tmp";
            const string& s1 = ss1.str();
            const string& s2 = ss2.str();
            remove (s1.c_str());
            remove (s2.c_str());
        }
    }

    // write time-domain signal
    const vec& ft_t1 = linspace(0.0,dt1*(nt1-1)/deom_fs2unit,nt1);
    const vec& ft_t2 = linspace(0.0,dt2*(nt2-1)/deom_fs2unit,nt2);
    const mat& ft_re = real(ft);
    const mat& ft_im = imag(ft);
    ft_t1.save("resp2nd.t1",raw_ascii);
    ft_t2.save("resp2nd.t2",raw_ascii);
    ft_re.save("resp2nd_re.t",raw_ascii);
    ft_im.save("resp2nd_im.t",raw_ascii);
    
    // 2D FFT
    ft = real(ft)*deom_c1;
    ft.row(0) *= 0.5;
    ft.col(0) *= 0.5;
    const cx_mat& fw2 = imag(ifft(ft))*nt2*dt2*deom_c1;
    const cx_mat& fw1 = ifft(fw2.st())*nt1*dt1;
    const cx_mat& fw = fw1.st();
    // const cx_mat& fw = ifft2(ft)*nt1*nt2*dt1*dt2;
    
    // write freq-domain signal
    const vec& fw_w1 = linspace(0.0,dw1*(nt1/2-1)/deom_cm2unit,nt1/2);
    const vec& fw_w2 = linspace(0.0,dw2*(nt2/2-1)/deom_cm2unit,nt2/2);
    const mat& fw_re = real(fw.submat(0,0,nt2/2-1,nt1/2-1));
    const mat& fw_im = imag(fw.submat(0,0,nt2/2-1,nt1/2-1));
    fw_w1.save("resp2nd.w1",raw_ascii);
    fw_w2.save("resp2nd.w2",raw_ascii);
    fw_re.save("resp2nd_re.w",raw_ascii);
    fw_im.save("resp2nd_im.w",raw_ascii);
}
}
