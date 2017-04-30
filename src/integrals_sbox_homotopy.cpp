/** \file integrals_sbox_homotopy.cpp

  Some driver/interface functions.    

  \todo Actually document this file.
*/

#include "integrals_sbox.h"
#include "integrals_sbox_formulation.h"

#include <memory>

extern "C" int luaopen_integrals_sbox_homotopy(lua_State*);

static int tau_gaussian_s(lua_State*);
static int tau_smoothbox_s(lua_State*);

static int tau_gaussian_symmetric_s(lua_State*);
static int tau_smoothbox_symmetric_s(lua_State*);

static int x_homotopy_s(lua_State*);
static int y_homotopy_s(lua_State*);

static int x_gaussian_b(lua_State*);
static int y_gaussian_b(lua_State*);
static int x_smoothbox_b(lua_State*);
static int y_smoothbox_b(lua_State*);

static int x_gaussian_symmetric_b(lua_State*);
static int x_smoothbox_symmetric_b(lua_State*);
static int y_gaussian_symmetric_b(lua_State*);
static int y_smoothbox_symmetric_b(lua_State*);

static inline int homotopy_s(lua_State*, cosine_or_sine);

static inline int tau_s(lua_State*,
                        const int,
                        tpgraphy_integrand_func); 

static inline int tau_symmetric_s(lua_State*,
                                  const int,
                                  tpgraphy_integrand_func);

static inline int gaussian_symmetric_b(lua_State*,
                                       const int,
                                       tpgraphy_theta_func,
                                       tpgraphy_integrand_func,
                                       removed_singularity_b_func);
   
extern "C" const char integrals_lualib_name[]   = "integrals_sbox";
extern "C" const char integrals_lualib_regkey[] = "integrals_sbox_rk";
extern "C" const luaL_Reg integralsR[] = {
    {"tau_gaussian_s", tau_gaussian_s},
    {"tau_smoothbox_s", tau_smoothbox_s},
    {"tau_gaussian_symmetric_s", tau_gaussian_symmetric_s},
    {"tau_smoothbox_symmetric_s", tau_smoothbox_symmetric_s},
    {"x_s", x_homotopy_s},
    {"y_s", y_homotopy_s},
    {"x_gaussian_b", x_gaussian_b},
    {"y_gaussian_b", y_gaussian_b},
    {"x_smoothbox_b", x_smoothbox_b},
    {"y_smoothbox_b", y_smoothbox_b},
    {"x_gaussian_symmetric_b", x_gaussian_symmetric_b},
    {"y_gaussian_symmetric_b", y_gaussian_symmetric_b},
    {"x_smoothbox_symmetric_b", x_smoothbox_symmetric_b},
    {"y_smoothbox_symmetric_b", y_smoothbox_symmetric_b},
    {"phi", pwlinear_cubic_homotopy_phi_luaf},
    {"dphi", pwlinear_cubic_homotopy_dphi_luaf},
    {NULL, NULL}
};

int luaopen_integrals_sbox_homotopy(lua_State *L)
{
    /** \todo May want to clarify the naming of the library. */
    return luaopen_integrals_sbox(L);
}

/* Is static ok here? */
static inline int homotopy_s(lua_State* L, cosine_or_sine cos_sin) {

    static const int NARG = 7;

    const int m = veclua_veclength(L, -NARG+0);
    const int n = veclua_veclength(L, -NARG+2);
    const int n_sub = veclua_veclength(L, -NARG+3);

    const std::unique_ptr<real[]> heap_v(new real[n+(n-1)+(2*m)+(3*n_sub)]);

    const real_vector theta_mid(    m,             &(heap_v[0]), L, -NARG+0);
    const real_vector   tau_mid(    m,             &(heap_v[m]), L, -NARG+1);
    const real_vector      beta(    n,           &(heap_v[2*m]), L, -NARG+2);
    const real_vector  beta_sub(n_sub,       &(heap_v[(2*m)+n]), L, -NARG+3);
    const real_vector   phi_sub(n_sub, &(heap_v[(2*m)+n+n_sub]), L, -NARG+4);
    const real_vector  dphi_sub(n_sub,
                                &(heap_v[(2*m)+n+(2*n_sub)]),
                                L,
                                -NARG+5); 
    /* output vector */
    real_vector z_s(n-1, &(heap_v[(2*m)+n+(3*n_sub)]));
  
    /* homotopy parameter */
    const real s = lua_tonumber(L, -NARG+6);

    lua_pop(L, NARG);

    real z0;
    z_s_worker(                    theta_mid,  tau_mid,
                                        beta, beta_sub,
                                     phi_sub, dphi_sub,
                                           s,
                pwlinear_cubic_homotopy_dphi,  cos_sin,                           
                                          z0,     z_s);

    veclua_pushtable(L, z_s);

    lua_pushnumber(L, z0);

    return 2;

}

/** Entry point to compute surface \f$x_s(\beta_{j+1/2})\f$ using \f$\tau\f$ and \f$\theta\f$.
  
   Arguments are passed by the lua stack
 - theta_mid \f$(n-1)\f$-vector; the angle of the free surface at \f$\beta_{j+1/2}\f$.
 - tau_mid \f$(n-1)\f$-vector; the log of the speed on the free surface at \f$\beta_{j+1/2}\f$.
 - beta \f$(n)\f$-vector; the grid points \f$\beta_j\f$.
 - beta_sub \f$(m)\f$-vector; \f$\beta\f$ at \f$\phi^{[j]}\f$.
 - phi_sub \f$(m)\f$-vector; 'clustering' locations in computational grid; \f$\phi^{[j]}\f$.
 - dphi_sub \f$(m)\f$-vector; value of \f$\phi'\f$ at \f$\phi^{[j]}\f$.
 - s scalar; optional homotopy parameter.
 \param L Lua stack pointer.
*/
int x_homotopy_s(lua_State* L) {
 
    return homotopy_s(L, USE_COSINE);

}

/** Entry point to compute surface \f$y_s(\beta_{j+1/2})\f$ using \f$\theta\f$ and \f$\tau\f$.
  
   Arguments are passed by the lua stack
 - theta_mid \f$(n-1)\f$-vector; the angle of the free surface at \f$\beta_{j+1/2}\f$.
 - tau_mid \f$(n-1)\f$-vector; the log of the speed on the free surface at \f$\beta_{j+1/2}\f$.
 - beta \f$(n)\f$-vector; the grid points \f$\beta_j\f$.
 - beta_sub \f$(m)\f$-vector; \f$\beta\f$ at \f$\phi^{[j]}\f$.
 - phi_sub \f$(m)\f$-vector; 'clustering' locations in computational grid; \f$\phi^{[j]}\f$.
 - dphi_sub \f$(m)\f$-vector; value of \f$\phi'\f$ at \f$\phi^{[j]}\f$.
 - s scalar; optional homotopy parameter.
 \param L Lua stack pointer.
*/
int y_homotopy_s(lua_State* L) {

    return homotopy_s(L, USE_SINE);

}


static inline int tau_s(lua_State*                        L,
                        const int                   n_t_arg,
                        tpgraphy_integrand_func integrand_b) {

    const int NARG = 11 + n_t_arg;

    const int m     = veclua_veclength(L, -NARG+0);
    const int n_sub = veclua_veclength(L, -NARG+1);
    const int n     = veclua_veclength(L, -NARG+3);

    const std::unique_ptr<real[]> heap_v(
                                      new real[m+n+(n-1)+(3*n_sub)+n_t_arg]
                                  );

    /* free surface grid arguments */
    const real_vector  theta_s(    m,             &(heap_v[0]), L, -NARG+0);
    const real_vector  phi_sub(n_sub,             &(heap_v[m]), L, -NARG+1);
    const real_vector dphi_sub(n_sub,           &(heap_v[m+n]), L, -NARG+2);
    const real_vector beta_sub(n_sub,     &(heap_v[m+n+n_sub]), L, -NARG+4);
    const real_vector     beta(    n, &(heap_v[m+n+(2*n_sub)]), L, -NARG+3);
    /* output vector */
    real_vector tau_s(n-1, &(heap_v[m+n+(3*n_sub)]));

    /* topography parameters */
    real_vector topography_params(4, &(heap_v[m+n+(n-1)+(3*n_sub)]));
    for(int j = 1; j <= n_t_arg; j++) {
        topography_params(j) = lua_tonumber(L, -NARG+4+j);
    }


    /* far field parameters */
    const real D_ffds      = lua_tonumber(L, -NARG+n_t_arg+5);
    const real D_ffus      = lua_tonumber(L, -NARG+n_t_arg+6);
    const real lambda_ffds = lua_tonumber(L, -NARG+n_t_arg+7);
    const real lambda_ffus = lua_tonumber(L, -NARG+n_t_arg+8);
    const real gamma       = lua_tonumber(L, -NARG+n_t_arg+9);

    /* homotopy parameter */
    const real s = lua_tonumber(L, -NARG+n_t_arg+10);

    lua_pop(L, NARG);

    ff_params ffus_params = {lambda_ffus,
                             gamma,
                             D_ffus,
                             phi_sub(1)};
    ff_params ffds_params = {lambda_ffds,
                             gamma,
                             D_ffds,
                             phi_sub(phi_sub.length())};

    lua_pushstring(L, integrals_lualib_regkey);
    lua_gettable(L, LUA_REGISTRYINDEX);
    gsl_integration_workspace **ws =
        (gsl_integration_workspace**)lua_touserdata(L, -1);

    method_s_funcs homotopy_method = {pwlinear_cubic_homotopy_phi,
                                      pwlinear_cubic_homotopy_dphi,
                                      compute_convolution_s,
                                      removed_singularity_general_s,
                                      far_upstream_integrand_s,
                                      far_downstream_integrand_s,
                                      integrand_b};

    tau_s_worker(          theta_s,
                           phi_sub,    dphi_sub,
                              beta,    beta_sub,
                 topography_params,
                       ffus_params, ffds_params,
                                 s,          ws,
                   homotopy_method,
                             tau_s);

    veclua_pushtable(L, tau_s);

    return 1;

}


int tau_gaussian_s(lua_State* L) {

    /* topography parameters
     * A - the amplitude
     * B - the coefficient of x in the exponential
     * L - the actual physical distance between inflection points */

    return tau_s(L, 3, gaussian_b);

}


int tau_smoothbox_s(lua_State* L) {

    /* topography parameters
     * A - the amplitude (depth)
     * B - rate (of decay) in the exponential terms
     * l - the actual physical distance between inflection points
     * phi_c - location of inflection points in phi */

    return tau_s(L, 4, smoothbox_b);

}


/** \todo documentation */
static inline int tau_symmetric_s(lua_State*                        L,
                                  const int                   n_t_arg,
                                  tpgraphy_integrand_func integrand_b) {

    const int NARG = 9 + n_t_arg;

    const int m     = veclua_veclength(L, -NARG+0);
    const int n_sub = veclua_veclength(L, -NARG+1);
    const int n     = veclua_veclength(L, -NARG+3);

    const std::unique_ptr<real[]> heap_v(
                                      new real[m+n+(n-1)+(3*n_sub)+n_t_arg]
                                  );

    /* free surface grid arguments */
    const real_vector  theta_s(    m,             &(heap_v[0]), L, -NARG+0);
    const real_vector  phi_sub(n_sub,             &(heap_v[m]), L, -NARG+1);
    const real_vector dphi_sub(n_sub,           &(heap_v[m+n]), L, -NARG+2);
    const real_vector beta_sub(n_sub,     &(heap_v[m+n+n_sub]), L, -NARG+4);
    const real_vector     beta(    n, &(heap_v[m+n+(2*n_sub)]), L, -NARG+3);
    /* output vector */
    real_vector tau_s(n-1, &(heap_v[m+n+(3*n_sub)]));

    /* topography parameters */
    real_vector topography_params(4, &(heap_v[m+n+(n-1)+(3*n_sub)]));
    for(int j = 1; j <= n_t_arg; j++) {
        topography_params(j) = lua_tonumber(L, -NARG+4+j);
    }

    /* far field parameters */
    const real D      = lua_tonumber(L, -NARG+n_t_arg+5);
    const real lambda = lua_tonumber(L, -NARG+n_t_arg+6);
    const real gamma  = lua_tonumber(L, -NARG+n_t_arg+7);

    /* homotopy parameter */
    const real s = lua_tonumber(L, -NARG+n_t_arg+8);

    lua_pop(L, NARG);

    ff_params ffus_params = {lambda,
                             gamma,
                             D,
                             -phi_sub(phi_sub.length())};
    ff_params ffds_params = {lambda,
                             gamma,
                             D,
                             phi_sub(phi_sub.length())};

    lua_pushstring(L, integrals_lualib_regkey);
    lua_gettable(L, LUA_REGISTRYINDEX);
    gsl_integration_workspace **ws =
        (gsl_integration_workspace**)lua_touserdata(L, -1);

    method_s_funcs homotopy_method = {pwlinear_cubic_homotopy_phi,
                                      pwlinear_cubic_homotopy_dphi,
                                      compute_convolution_symmetric_s,
                                      removed_singularity_symmetric_s,
                                      far_upstream_integrand_s,
                                      far_downstream_integrand_s,
                                      integrand_b};

    tau_s_worker(          theta_s,
                           phi_sub,    dphi_sub, 
                              beta,    beta_sub,
                 topography_params,
                       ffus_params, ffds_params,
                                 s,          ws,
                   homotopy_method,
                             tau_s);

    veclua_pushtable(L, tau_s);

    return 1;

}


int tau_gaussian_symmetric_s(lua_State* L) {

    /* topography parameters
     * A - the amplitude
     * B - the coefficient of x in the exponential
     * L - the actual physical distance between inflection points */

    return tau_symmetric_s(L, 3, gaussian_b);

}


int tau_smoothbox_symmetric_s(lua_State* L) {

    /* topography parameters
     * A - the amplitude (depth)
     * B - rate (of decay) in the exponential terms
     * l - the actual physical distance between inflection points
     * phi_c - location of inflection points in phi */

    return tau_symmetric_s(L, 4, smoothbox_b);

}


static inline int z_b(lua_State* L,
                      const int n_t_arg,
                      tpgraphy_theta_func theta_f,
                      tpgraphy_integrand_func integrand_b,
                      removed_singularity_b_func singularity_f,
                      cosine_or_sine cos_sin) {

    const int NARG = 12 + n_t_arg;

    const int m     = veclua_veclength(L, -NARG+0);
    const int n     = veclua_veclength(L, -NARG+1);
    const int n_sub = veclua_veclength(L, -NARG+3);

    const std::unique_ptr<real[]> heap_v(
                                      new real[(2*m)+(2*n)+(3*n_sub)+n_t_arg]
                                  );

    /* free surface grid arguments */
    const real_vector      phi(    m,                 &(heap_v[0]), L, -NARG+0);
    const real_vector  theta_s(    n,                 &(heap_v[m]), L, -NARG+1);
    const real_vector   beta_s(    n,               &(heap_v[m+n]), L, -NARG+2);
    const real_vector  phi_sub(n_sub,           &(heap_v[m+(2*n)]), L, -NARG+3);
    const real_vector dphi_sub(n_sub,     &(heap_v[m+(2*n)+n_sub]), L, -NARG+4);
    const real_vector beta_sub(n_sub, &(heap_v[m+(2*n)+(2*n_sub)]), L, -NARG+5);
    /* output vector */
    real_vector z_b(m, &(heap_v[m+(2*n)+(3*n_sub)]));

    /* topography parameters */
    real_vector topography_params(4, &(heap_v[(2*m)+(2*n)+(3*n_sub)]));
    for(int j = 1; j <= n_t_arg; j++) {
        topography_params(j) = lua_tonumber(L, -NARG+5+j);
    }

    /* far field parameters */
    const real D_ffus      = lua_tonumber(L, -NARG+n_t_arg+6);
    const real D_ffds      = lua_tonumber(L, -NARG+n_t_arg+7);
    const real lambda_ffus = lua_tonumber(L, -NARG+n_t_arg+8);
    const real lambda_ffds = lua_tonumber(L, -NARG+n_t_arg+9);
    const real gamma  = lua_tonumber(L, -NARG+n_t_arg+10);

    /* homotopy parameter */
    const real s = lua_tonumber(L, -NARG+n_t_arg+11);


    ff_params ffus_params = {lambda_ffus,
                             gamma,
                             D_ffus,
                             phi_sub(1)};
    ff_params ffds_params = {lambda_ffds,
                             gamma,
                             D_ffds,
                             phi_sub(phi_sub.length())};

    lua_pushstring(L, integrals_lualib_regkey);
    lua_gettable(L, LUA_REGISTRYINDEX);
    gsl_integration_workspace **ws = 
        (gsl_integration_workspace**)lua_touserdata(L, -1);

    method_b_funcs homotopy_method = {pwlinear_cubic_homotopy_phi,
                                      pwlinear_cubic_homotopy_dphi,
                                      theta_f,
                                      topography_compute_convolution_s,
                                      integrand_b,
                                      singularity_f,
                                      topography_far_upstream_s,
                                      topography_far_downstream_s};

    z_b_worker(              phi,
                         theta_s,      beta_s,
                        beta_sub,     phi_sub, dphi_sub,
                               s,          ws,
                 homotopy_method,
               topography_params,
                     ffus_params, ffds_params,
                         cos_sin,
                             z_b);

    veclua_pushtable(L, z_b);

    return 1;

}


int x_gaussian_b(lua_State* L) {

    return z_b(L,
               3,
               gaussian_theta,
               topography_gaussian_integrand_b,
               topography_gaussian_removed_singularity_b,
               USE_COSINE);

}


int y_gaussian_b(lua_State* L) {

    return z_b(L,
               3,
               gaussian_theta,
               topography_gaussian_integrand_b,
               topography_gaussian_removed_singularity_b,
               USE_SINE);

}


int x_smoothbox_b(lua_State* L) {

    return z_b(L,
               4,
               smoothbox_theta,
               topography_smoothbox_integrand_b,
               topography_smoothbox_removed_singularity_b,
               USE_COSINE);

}


int y_smoothbox_b(lua_State* L) {

    return z_b(L,
               4,
               smoothbox_theta,
               topography_smoothbox_integrand_b,
               topography_smoothbox_removed_singularity_b,
               USE_SINE);

}


inline int z_symmetric_b(lua_State* L,
                         const int n_t_arg,
                         tpgraphy_theta_func theta_f,
                         tpgraphy_integrand_func integrand_b,
                         removed_singularity_b_func singularity_f,
                         cosine_or_sine cos_sin) {

    const int NARG = 10 + n_t_arg;

    const int m     = veclua_veclength(L, -NARG+0);
    const int n     = veclua_veclength(L, -NARG+1);
    const int n_sub = veclua_veclength(L, -NARG+3);

    const std::unique_ptr<real[]> heap_v(
                                      new real[(2*m)+(2*n)+(3*n_sub)+n_t_arg]
                                  );

    /* free surface grid arguments */
    const real_vector      phi(    m,                 &(heap_v[0]), L, -NARG+0);
    const real_vector  theta_s(    n,                 &(heap_v[m]), L, -NARG+1);
    const real_vector   beta_s(    n,               &(heap_v[m+n]), L, -NARG+2);
    const real_vector  phi_sub(n_sub,           &(heap_v[m+(2*n)]), L, -NARG+3);
    const real_vector dphi_sub(n_sub,     &(heap_v[m+(2*n)+n_sub]), L, -NARG+4);
    const real_vector beta_sub(n_sub, &(heap_v[m+(2*n)+(2*n_sub)]), L, -NARG+5);
    /* output vector */
    real_vector z_b(m, &(heap_v[m+(2*n)+(3*n_sub)]));

    /* topography parameters */
    real_vector topography_params(4, &(heap_v[(2*m)+(2*n)+(3*n_sub)]));
    for(int j = 1; j <= n_t_arg; j++) {
        topography_params(j) = lua_tonumber(L, -NARG+5+j);
    }

    /* far field parameters */
    const real D      = lua_tonumber(L, -NARG+n_t_arg+6);
    const real lambda = lua_tonumber(L, -NARG+n_t_arg+7);
    const real gamma  = lua_tonumber(L, -NARG+n_t_arg+8);

    /* homotopy parameter */
    const real s = lua_tonumber(L, -NARG+n_t_arg+9);


    ff_params ffus_params = {lambda, gamma, D, -phi_sub(phi_sub.length())};
    ff_params ffds_params = {lambda, gamma, D,  phi_sub(phi_sub.length())};

    lua_pushstring(L, integrals_lualib_regkey);
    lua_gettable(L, LUA_REGISTRYINDEX);
    gsl_integration_workspace **ws = 
        (gsl_integration_workspace**)lua_touserdata(L, -1);

    method_b_funcs homotopy_method = {pwlinear_cubic_homotopy_phi,
                                      pwlinear_cubic_homotopy_dphi,
                                      theta_f,
                                      topography_compute_symmetric_convolution_s,
                                      integrand_b,
                                      singularity_f,
                                      topography_far_upstream_s,
                                      topography_far_downstream_s};

    z_b_worker(              phi,
                         theta_s,      beta_s,
                        beta_sub,     phi_sub, dphi_sub,
                               s,          ws,
                 homotopy_method,
               topography_params,
                     ffus_params, ffds_params,
                         cos_sin,
                             z_b);

    veclua_pushtable(L, z_b);

    return 1;

}


int x_gaussian_symmetric_b(lua_State* L) {

    return z_symmetric_b(L,
                         3,
                         gaussian_theta,
                         topography_gaussian_integrand_b,
                         topography_gaussian_removed_singularity_b,
                         USE_COSINE);

}


int y_gaussian_symmetric_b(lua_State* L) {

    return z_symmetric_b(L,
                         3,
                         gaussian_theta,
                         topography_gaussian_integrand_b,
                         topography_gaussian_removed_singularity_b,
                         USE_SINE);

}


int x_smoothbox_symmetric_b(lua_State* L) {

    return z_symmetric_b(L,
                         4,
                         smoothbox_theta,
                         topography_smoothbox_integrand_b,
                         topography_smoothbox_removed_singularity_b,
                         USE_COSINE);

}


int y_smoothbox_symmetric_b(lua_State* L) {

    return z_symmetric_b(L,
                         4,
                         smoothbox_theta,
                         topography_smoothbox_integrand_b,
                         topography_smoothbox_removed_singularity_b,
                         USE_SINE);

}

