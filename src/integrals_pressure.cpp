/** \file integrals_pressure.cpp

  C implementation of integral formulas.

*/

#include "integrals_pressure.h"

#include <memory>

static int integrals_pressure_clean_up(lua_State *);

extern char integrals_lualib_name[];
extern char integrals_lualib_regkey[];

extern luaL_Reg integralsR[];

int luaopen_integrals_pressure(lua_State *L) {
    /* Create userdata with the gsl_integration_workspace */
    lua_newtable(L);

    const luaL_Reg* fs = integralsR;

    while(fs->name) {
        lua_pushstring(L, fs->name);
        lua_pushcfunction(L, fs->func);
        lua_settable(L, -3);
        ++fs;
    }
 
    /** \todo Fix which name to use for the gsl_integration_workspace registry
     * 'key' */
    lua_pushstring(L, integrals_lualib_regkey);
    // Make new userdata at the top of the stack
    gsl_integration_workspace **ws =
        (gsl_integration_workspace**)lua_newuserdata(
                                         L,
                                         sizeof(gsl_integration_workspace **)
                                     );
    *ws = gsl_integration_workspace_alloc(1000);

    // Make a new metatable
    lua_newtable(L);
    // Garbage collection method
    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, integrals_pressure_clean_up);
    lua_settable(L, -3);
    // Set metatable
    lua_setmetatable(L, -2);
  
    lua_settable(L, LUA_REGISTRYINDEX);
  
    return 1;

}

/* Clean up function for the gsl_integration_workspace */
int integrals_pressure_clean_up(lua_State* L) {

    lua_pushstring(L, integrals_lualib_regkey);
    lua_gettable(L, LUA_REGISTRYINDEX);
    gsl_integration_workspace **ws = 
        (gsl_integration_workspace**)lua_touserdata(L, -1);
  
    gsl_integration_workspace_free(*ws);

    return 0;

}


void u_s_worker(const real_vector&         v_s,
                const real_vector&     phi_sub,
                const real_vector&    dphi_sub,
                const real_vector&        beta,
                const real_vector &   beta_sub,
                const ff_params &  ffus_params,
                const ff_params &  ffds_params,
                const real                   s,
                gsl_integration_workspace** ws,
                const method_funcs &    method,
                real_vector &              u_s) {

    const unsigned int n = beta.length();

    /* Heap for vector needed 
     * name           size
     * phi             (n)
     * dphi            (n)
     * integrand       (n)
     * beta_mid      (n-1)
     * phi_mid       (n-1)
     * singularity   (n-1)
     * ffds_integral (n-1)
     * ffus_integral (n-1) */
    const std::unique_ptr<real[]> heap_v(new real[(3*n)+(5*(n-1))]);
    
    /* grid parameters */
    real_vector phi(n, &(heap_v[0]));
    real_vector dphi(n, &(heap_v[n]));
    
    real_vector beta_mid(n-1, &(heap_v[2*n]));
    real_vector phi_mid(n-1, &(heap_v[(2*n)+(n-1)]));

    /* convolution integrand */
    real_vector integrand(n, &(heap_v[2*n+(2*(n-1))]));
    
    /* removed singularity */
    real_vector singularity(n-1, &(heap_v[(3*n)+(2*(n-1))]));
    
    /* far field integrals */
    real_vector ffds_integral(n-1, &(heap_v[(3*n)+(3*(n-1))]));
    real_vector ffus_integral(n-1, &(heap_v[(3*n)+(4*(n-1))]));
    
    for(unsigned int i = 1; i <= n-1; i++) {
        beta_mid(i) = (beta(i) + beta(i+1)) / 2.0;
    }

    method.transform_phi(beta, beta_sub, phi_sub, dphi_sub, phi, s);
    method.transform_dphi(beta, beta_sub, phi_sub, dphi_sub, dphi, s);
  
    method.transform_phi(beta_mid, beta_sub, phi_sub, dphi_sub, phi_mid, s);
  
    method.removed_singularity(v_s, phi, phi_mid, singularity);

    method.linear_downstream_int(phi_mid,
                                 ffds_params.phi_match,
                                 ffds_params.D,
                                 ffds_params.lambda,
                                 ffds_integral);

    method.linear_upstream_int(phi_mid,
                               ffus_params.phi_match,
                               ffus_params.D,
                               ffus_params.lambda,
                               ffus_integral);

    for(unsigned int i = 1; i <= n-1; i++) {
        /* From the formulation, begin with u = 1.0 */
        u_s(i) = 1.0;

        method.convolution_terms(i, v_s, phi, dphi, phi_mid, integrand);

        for(unsigned int j = 1; j <= n; j++) {
            const real delta = (j == 1 ? beta(2) - beta(1) :
                                    (j == n ? beta(n) - beta(n-1) :
                                         beta(j+1) - beta(j-1)));
            u_s(i) += delta * integrand(j) / 2.0;
        }

        u_s(i) += singularity(i);
        u_s(i) += ffds_integral(i);
        u_s(i) += ffus_integral(i);
    }

}


void z_s_worker(const real_vector &     u_mid,
                const real_vector &     v_mid,
                const real_vector &      beta,
                const real_vector &  beta_sub,
                const real_vector &   phi_sub,
                const real_vector &  dphi_sub,
                const real                  s,
                transform_func transform_dphi,
                const use_u_or_v       u_or_v,
                real &                     z0,
                real_vector &             z_s) {

    const int n = beta.length();

    const std::unique_ptr<real[]> heap_v(new real[2*(n-1)]);

    real_vector beta_mid(n-1, &(heap_v[0]));
    real_vector dphi_mid(n-1, &(heap_v[n-1]));

    for(unsigned int i = 1; i <= n-1; i++) {
        beta_mid(i) = (beta(i) + beta(i+1)) / 2.0;
    }

    transform_dphi(beta_mid, beta_sub, phi_sub, dphi_sub, dphi_mid, s);

    real M_j = 0.0;

    for(unsigned int j = n-1; j >= 1; j--) {
        const real f_j = -dphi_mid(j) *
                             (u_or_v == USE_V ? v_mid(j) : u_mid(j)) /
                             (u_mid(j) * u_mid(j) + v_mid(j) * v_mid(j));
        M_j += f_j * (beta(j+1) - beta(j));

        z_s(j) = M_j - f_j * (beta(j+1) - beta(j)) / 2.0;
    }

    z0 = M_j;

}


