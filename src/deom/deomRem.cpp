/**
 * git clone https://github.com/hou-dao/raman.git
 * ---
 * Writtened by Houdao Zhang 
 * mailto: houdao@connect.ust.hk
 */
#include <cmath>
#include "deom.hpp"

void deom::rem_sch (cx_cube& dtotal, const cx_cube& total, const double t) {

    const int nsav = nddo;
    dtotal.slices(0,nddo-1).zeros();

    for (int iado=0; iado<nsav; ++iado) {
        const cx_mat& ado = total.slice(iado);
        if (iado==0 || is_valid (ado)) {
            string key = keys(iado);
            const int tier = maps[key].tier;
            const cx_double gams = maps[key].gams;

            dtotal.slice(iado) += -deom_ci*(ham1*ado-ado*ham1)-gams*ado;
            for (int m=0; m<nmod; ++m) {
                qddo.slice(m) = qmd1.slice(m)*ado;
                ddoq.slice(m) = ado*qmd1.slice(m);
                dtotal.slice(iado) -= delt_res(m)*(
                        qmd1.slice(m)*(qddo.slice(m)-ddoq.slice(m))
                       -(qddo.slice(m)-ddoq.slice(m))*qmd1.slice(m));
            }

            if (tier < lmax) {
                for (int mp=0; mp<nind; ++mp) {
                    const int m = mp/nper;
                    const int p = mp%nper;
                    const int n = key[mp]-'0';
                    const cx_double sn = -deom_ci*sqrt((n+1)/coef_abs(m,p));
                    const cx_double cl = sn*coef_lft(m,p);
                    const cx_double cr = sn*coef_rht(m,p);
                    key[mp] = static_cast<char>(key[mp]+1);
                    const auto &xp = maps.emplace(key,hnod(nddo,tier+1,gams+expn(mp)));
                    if (!xp.second) {
                        int loc = xp.first->second.rank;
                        dtotal.slice(loc) += cl*qddo.slice(m)-cr*ddoq.slice(m);
                    } else {
                        keys(nddo) = key;
                        dtotal.slice(nddo) = cl*qddo.slice(m)-cr*ddoq.slice(m);
                        nddo += 1;
                    }
                    key[mp] = static_cast<char>(key[mp]-1);
                }
            }

            if (iado > 0) {
                for (int mp=0; mp<nind; ++mp) {
                    const int m = mp/nper;
                    const int p = mp%nper;
                    const int n = key[mp]-'0';
                    if (n > 0) {
                        key[mp] = static_cast<char>(key[mp]-1);
                        const cx_double sn = -deom_ci*sqrt(n*coef_abs(m,p));
                        const auto &xp = maps.emplace(key,hnod(nddo,tier-1,gams-expn(mp)));
                        if (!xp.second) {
                            int loc = xp.first->second.rank;
                            dtotal.slice(loc) += sn*(qddo.slice(m)-ddoq.slice(m)); 
                        } else {
                            keys(nddo) = key;
                            dtotal.slice(nddo) = sn*(qddo.slice(m)-ddoq.slice(m)); 
                            nddo += 1;
                        }
                        key[mp] = static_cast<char>(key[mp]+1);
                    }
                }
            }
        }
    }
}

void deom::rem_hei (cx_cube& dtotal, const cx_cube& total, const double t) {

    const int nsav = nddo;
    dtotal.slices(0,nddo-1).zeros();

    for (int iado=0; iado<nsav; ++iado) {
        const cx_mat& ado = total.slice(iado);
        if (iado==0 || is_valid (ado)) {
            string key = keys(iado);
            const int tier = maps[key].tier;
            const cx_double gams = maps[key].gams;

            dtotal.slice(iado) += -deom_ci*(ado*ham1-ham1*ado)-gams*ado;
            for (int m=0; m<nmod; ++m) {
                qddo.slice(m) = qmd1.slice(m)*ado;
                ddoq.slice(m) = ado*qmd1.slice(m);
                dtotal.slice(iado) -= delt_res(m)*(
                        qmd1.slice(m)*(qddo.slice(m)-ddoq.slice(m))
                       -(qddo.slice(m)-ddoq.slice(m))*qmd1.slice(m));
            }

            if (iado > 0) {
                for (int mp=0; mp<nind; ++mp) {
                    const int m = mp/nper;
                    const int p = mp%nper;
                    const int n = key[mp]-'0';
                    if (n > 0) {
                        const cx_double sn = -deom_ci*sqrt(n/coef_abs(m,p));
                        const cx_double cl = sn*coef_lft(m,p);
                        const cx_double cr = sn*coef_rht(m,p);
                        key[mp] = static_cast<char>(key[mp]-1);
                        const auto &xp = maps.emplace(key,hnod(nddo,tier-1,gams-expn(mp)));
                        if (!xp.second) {
                            int loc = xp.first->second.rank;
                            dtotal.slice(loc) += cl*ddoq.slice(m)-cr*qddo.slice(m);
                        } else {
                            keys(nddo) = key;
                            dtotal.slice(nddo) = cl*ddoq.slice(m)-cr*qddo.slice(m);
                            nddo += 1;
                        }
                        key[mp] = static_cast<char>(key[mp]+1);
                    }
                }
            }

            if (tier < lmax) {
                for (int mp=0; mp<nind; ++mp) {
                    const int m = mp/nper;
                    const int p = mp%nper;
                    const int n = key[mp]-'0';
                    key[mp] = static_cast<char>(key[mp]+1);
                    const cx_double sn = -deom_ci*sqrt((n+1)*coef_abs(m,p));
                    const auto &xp = maps.emplace(key,hnod(nddo,tier+1,gams+expn(mp)));
                    if (!xp.second) {
                        int loc = xp.first->second.rank;
                        dtotal.slice(loc) += sn*(ddoq.slice(m)-qddo.slice(m)); 
                    } else {
                        keys(nddo) = key;
                        dtotal.slice(nddo) = sn*(ddoq.slice(m)-qddo.slice(m)); 
                        nddo += 1;
                    }
                    key[mp] = static_cast<char>(key[mp]-1);
                }
            }
        }
    }
}
