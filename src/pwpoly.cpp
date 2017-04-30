/** \file pwpoly.cpp

    Functions for managing piecewise polynomial transformations of domain

*/

#include "pwpoly.h"
#include <memory>

static inline int pwlinear_cubic_homotopy_luaf(lua_State*,
                                               const transform_func);

void pwlinear_cubic_homotopy_phi(const real_vector&     beta,
                                 const real_vector& beta_sub,
                                 const real_vector&  phi_sub,
                                 const real_vector& dphi_sub,
                                 real_vector              &o,
                                 const real                s) {

    const unsigned int m = phi_sub.length();

    unsigned int j = 1;
    for(unsigned int i = 1; i <= beta.length(); i++) {

        if(beta(i) <= beta_sub(1)) {
            o(i) = phi_sub(1) + (beta(i) - beta_sub(1)) * dphi_sub(1);
        } else if(beta(i) > beta_sub(m)) {
            o(i) = phi_sub(m) + (beta(i) - beta_sub(m)) * dphi_sub(m);
        } else {
            while(beta_sub(j) < beta(i))
                j++;
      
            // beta_sub(j) is now right-hand end point, beta_sub(j-1) is
            // left-hand end point

            // Okay so the following terms help with the factors of the basis
            // polynomials
            const real dbeta = beta_sub(j) - beta_sub(j-1);

            const real bi_sub_bj0 = beta(i) - beta_sub(j-1);
            const real bi_sub_bj1 = beta(i) - beta_sub(j);

            const real zero1 = (3.0*beta_sub(j-1) - beta_sub(j)) / 2.0;
            const real zero2 = (3.0*beta_sub(j) - beta_sub(j-1)) / 2.0;

            // Calculate values of basis polynomials
            const real f1 = bi_sub_bj1 * bi_sub_bj1 *
                                (beta(i) - zero1) /
                                (dbeta*dbeta*(beta_sub(j-1) - zero1));
            const real f2 = bi_sub_bj0 * bi_sub_bj0 *
                                (beta(i) - zero2) /
                                (dbeta*dbeta*(beta_sub(j) - zero2));
            const real f3 = bi_sub_bj0 * bi_sub_bj1 * 
                                bi_sub_bj1 / (dbeta*dbeta);
            const real f4 = bi_sub_bj1 * bi_sub_bj0 * 
                                bi_sub_bj0 / (dbeta*dbeta);

            // Construct phi using homotopy between linear interpolation and the
            // cubic interpolation with fixed, continuous, derivatives
            o(i) = (1.0-s) *
                       (phi_sub(j)*bi_sub_bj0 - phi_sub(j-1)*bi_sub_bj1) /
                       dbeta;
            o(i) += s*(f1*phi_sub(j-1)  + f2*phi_sub(j) +
                           f3*dphi_sub(j-1) + f4*dphi_sub(j));
        }

    }

}

void pwlinear_cubic_homotopy_dphi(const real_vector&     beta,
                                  const real_vector& beta_sub,
                                  const real_vector&  phi_sub,
                                  const real_vector& dphi_sub,
                                  real_vector              &o,
                                  const real                s) {

    unsigned int m = phi_sub.length();

    unsigned int j = 1;
    for(unsigned int i = 1; i <= beta.length(); i++) {

        if(beta(i) <= beta_sub(1)) {
            o(i) = dphi_sub(1);
        } else if(beta(i) > beta_sub(m)) {
            o(i) = dphi_sub(m);
        } else {
            while(beta_sub(j) < beta(i))
                j++;

            // beta_sub(j) is now right-hand end point, beta_sub(j-1) is
            // left-hand end point

            // Okay so the following terms help with the factors of the basis
            // polynomials
            const real dbeta = beta_sub(j) - beta_sub(j-1);

            const real bi_sub_bj0 = beta(i) - beta_sub(j-1);
            const real bi_sub_bj1 = beta(i) - beta_sub(j);

            const real zero1 = (3.0*beta_sub(j-1) - beta_sub(j)) / 2.0;
            const real zero2 = (3.0*beta_sub(j) - beta_sub(j-1)) / 2.0;

            // Calculate values of basis polynomials
            const real df1 = bi_sub_bj1 * (2.0*(beta(i) - zero1) + bi_sub_bj1) /
                                 (dbeta*dbeta*(beta_sub(j-1) - zero1));
            const real df2 = bi_sub_bj0 * (2.0*(beta(i) - zero2) + bi_sub_bj0) /
                                 (dbeta*dbeta*(beta_sub(j) - zero2));
            const real df3 = bi_sub_bj1*(2.0*bi_sub_bj0 + bi_sub_bj1) /
                                 (dbeta*dbeta);
            const real df4 = bi_sub_bj0*(2.0*bi_sub_bj1 + bi_sub_bj0) /
                                 (dbeta*dbeta);

            // Construct phi using homotopy between linear interpolation and the
            // cubic interpolation with fixed, continuous, derivatives
            o(i)  = (1.0-s) * (phi_sub(j) - phi_sub(j-1)) / (dbeta);
            o(i) += s*(df1*phi_sub(j-1)  + df2*phi_sub(j) +
                           df3*dphi_sub(j-1) + df4*dphi_sub(j));
        }

    }

}

/* is static ok here? */
static inline int pwlinear_cubic_homotopy_luaf(lua_State* L,
                                               const transform_func f) {

    static const int NARG = 5;

    const int n     = veclua_veclength(L, -NARG+0);
    const int n_sub = veclua_veclength(L, -NARG+1);
    
    const std::unique_ptr<real[]> heap_v(new real[(2*n)+(3*n_sub)]);

    const real_vector     beta(    n,           &(heap_v[0]), L, -NARG+0);
    const real_vector beta_sub(n_sub,           &(heap_v[n]), L, -NARG+1);
    const real_vector  phi_sub(n_sub,     &(heap_v[n+n_sub]), L, -NARG+2);
    const real_vector dphi_sub(n_sub, &(heap_v[n+(2*n_sub)]), L, -NARG+3);
    /* output vector */
    real_vector o(n, &(heap_v[n+(3*n_sub)]));
//    veclua_tovector(L, -NARG+0, beta);
//    veclua_tovector(L, -NARG+1, beta_sub);
//    veclua_tovector(L, -NARG+2, phi_sub);
//    veclua_tovector(L, -NARG+3, dphi_sub);
  
    const real s = lua_tonumber(L, -NARG+4);

    lua_pop(L, NARG);

    f(beta, beta_sub, phi_sub, dphi_sub, o, s);

    veclua_pushtable(L, o);

    return 1;

}

int pwlinear_cubic_homotopy_phi_luaf(lua_State* L) {

    return pwlinear_cubic_homotopy_luaf(L, pwlinear_cubic_homotopy_phi);

}

int pwlinear_cubic_homotopy_dphi_luaf(lua_State* L) {

    return pwlinear_cubic_homotopy_luaf(L, pwlinear_cubic_homotopy_dphi);

}


