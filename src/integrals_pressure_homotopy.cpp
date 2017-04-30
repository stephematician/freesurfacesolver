/** \file integrals_pressure_homotopy.cpp
  
*/

#include "integrals_pressure.h"
#include "integrals_pressure_formulation.h"

#include <memory>

extern "C" int luaopen_integrals_pressure_homotopy(lua_State*);

/* General formulation */
static int x_homotopy_s(lua_State*);
static int y_homotopy_s(lua_State*);
static int u_homotopy_s(lua_State*);

/* Symmetric formulation */
static int u_homotopy_symmetric_s(lua_State*);

static inline int homotopy_s(lua_State*, const use_u_or_v);

extern "C" const char integrals_lualib_name[]   = "integrals_pressure";
extern "C" const char integrals_lualib_regkey[] = "integrals_pressure_rk";
extern "C" const luaL_Reg integralsR[] = {
    {"u_s", u_homotopy_s},
    {"x_s", x_homotopy_s},
    {"y_s", y_homotopy_s},
    {"u_symmetric_s", u_homotopy_symmetric_s},
    {"phi", pwlinear_cubic_homotopy_phi_luaf},
    {"dphi", pwlinear_cubic_homotopy_dphi_luaf},
    {NULL, NULL}
};

int luaopen_integrals_pressure_homotopy(lua_State *L) {

    /** \todo May want to clarify the naming of the library. */
    return luaopen_integrals_pressure(L);

}

/** Entry point to compute surface \f$x_s(\beta_{j+1/2})\f$ using \f$u\f$ and
    \f$v\f$.

    Arguments are passed by the Lua stack
    - u_mid \f$(n-1)\f$-vector; the horizontal component of velocity on the free
            surface at \f$\beta_{j+1/2}\f$.
    - v_mid \f$(n-1)\f$-vector; the vertical component of velocity on the free
            surface at \f$\beta_{j+1/2}\f$.
    - beta \f$(n)\f$-vector; the grid points \f$\beta_j\f$.
    - beta_sub \f$(m)\f$-vector; \f$\beta\f$ at \f$\phi^{[j]}\f$.
    - phi_sub \f$(m)\f$-vector; 'clustering' locations in computational grid;
              \f$\phi^{[j]}\f$.
    - dphi_sub \f$(m)\f$-vector; value of \f$\phi'\f$ at \f$\phi^{[j]}\f$.
    - s scalar; optional homotopy parameter.
    
    \param L Lua stack pointer.

*/
int x_homotopy_s(lua_State* L) {
    
    return homotopy_s(L, USE_U);

}


/** Entry point to compute surface \f$y_s(\beta_{j+1/2})\f$ using \f$u\f$ and
    \f$v\f$.
  
    Arguments are passed by the lua stack
    - u_mid \f$(n-1)\f$-vector; the horizontal component of velocity on the free
            surface at \f$\beta_{j+1/2}\f$.
    - v_mid \f$(n-1)\f$-vector; the vertical component of velocity on the free
            surface at \f$\beta_{j+1/2}\f$.
    - beta \f$(n)\f$-vector; the grid points \f$\beta_j\f$.
    - beta_sub \f$(m)\f$-vector; \f$\beta\f$ at \f$\phi^{[j]}\f$.
    - phi_sub \f$(m)\f$-vector; 'clustering' locations in computational grid;
              \f$\phi^{[j]}\f$.
    - dphi_sub \f$(m)\f$-vector; value of \f$\phi'\f$ at \f$\phi^{[j]}\f$.
    - s scalar; optional homotopy parameter.

    \param L Lua stack pointer.

*/
int y_homotopy_s(lua_State* L) {
    
    return homotopy_s(L, USE_V);

}

static inline int homotopy_s(lua_State* L, const use_u_or_v u_or_v) {

    static const int NARG = 7;

    const int m = veclua_veclength(L, -NARG+0);
    const int n = veclua_veclength(L, -NARG+2);
    const int n_sub = veclua_veclength(L, -NARG+3);

    std::unique_ptr<real[]> heap_v(new real[(3*m)+(2*n)+(3*n_sub)]);

    const real_vector    u_mid(    m,                 &(heap_v[0]), L, -NARG+0);
    const real_vector    v_mid(    m,                 &(heap_v[m]), L, -NARG+1);
    const real_vector     beta(    n,               &(heap_v[2*m]), L, -NARG+2);
    const real_vector beta_sub(n_sub,           &(heap_v[(2*m)+n]), L, -NARG+3);
    const real_vector  phi_sub(n_sub,     &(heap_v[(2*m)+n+n_sub]), L, -NARG+4);
    const real_vector dphi_sub(n_sub, &(heap_v[(2*m)+n+(2*n_sub)]), L, -NARG+5);
    
    real s = lua_tonumber(L, -NARG+6);

    lua_pop(L, NARG);
    
    real_vector      z_s(  n-1, &(heap_v[(2*m)+n+(3*n_sub)]));
    real z0;

  
    z_s_worker(u_mid, v_mid,
               beta,
               beta_sub, phi_sub, dphi_sub,
               s,
               pwlinear_cubic_homotopy_dphi,
               u_or_v,
               z0,
               z_s);

    veclua_pushtable(L, z_s);

    lua_pushnumber(L, z0);

    return 2;

}


/** Entry point to compute surface \f$u_s(\beta_{j+1/2})\f$ using \f$v\f$.
  
    Arguments are passed by the lua stack
    - v_s \f$(n-1)\f$-vector; the horizontal component of velocity on the free
          surface at \f$\beta_{j+1/2}\f$.
    - phi_sub \f$(m)\f$-vector; 'clustering' locations in computational grid;
              \f$\phi^{[j]}\f$.
    - dphi_sub \f$(m)\f$-vector; value of \f$\phi'\f$ at \f$\phi^{[j]}\f$.
    - beta \f$(n)\f$-vector; the grid points \f$\beta_j\f$.
    - beta_sub \f$(m)\f$-vector; \f$\beta\f$ at \f$\phi^{[j]}\f$.
    - D_ffds scalar; coefficent of linearised far field downstream velocity.
    - D_ffus scalar; coefficent of linearised far field upstream velocity.
    - lambda_ffds scalar; exponent of linearised far field downstream velocity.
    - lambda_ffus scalar; exponent of linearised far field upstream velocity.
    - s scalar; optional homotopy parameter.

 \param L Lua stack pointer.
 
*/
int u_homotopy_s(lua_State* L) {

    static const int NARG = 10;

    const int n = veclua_veclength(L, -NARG+0);
    const int n_sub = veclua_veclength(L, -NARG+1);

    const std::unique_ptr<real[]> heap_v(new real [(2*n)+(n-1)+(3*n_sub)]);

    /* free surface variables */
    const real_vector v_s(n, &(heap_v[0]), L, -NARG+0);
    /* grid variables */
    const real_vector     beta(    n,               &(heap_v[n]), L, -NARG+3);
    const real_vector  phi_sub(n_sub,             &(heap_v[2*n]), L, -NARG+1);
    const real_vector dphi_sub(n_sub,     &(heap_v[(2*n)+n_sub]), L, -NARG+2);
    const real_vector beta_sub(n_sub, &(heap_v[(2*n)+(2*n_sub)]), L, -NARG+4);

    /* far field parameters */
    const real D_ffds      = lua_tonumber(L, -NARG+5);
    const real D_ffus      = lua_tonumber(L, -NARG+6);
    const real lambda_ffds = lua_tonumber(L, -NARG+7);
    const real lambda_ffus = lua_tonumber(L, -NARG+8);

    /* homotopy parameter */
    const real s = lua_tonumber(L, -NARG+9);

    lua_pop(L, NARG);

    real_vector u_s(n, &(heap_v[(2*n)+(3*n_sub)]));

    ff_params ffus_params = {lambda_ffus, 0, D_ffus, phi_sub(1)};
    ff_params ffds_params = {lambda_ffds, 0, D_ffds, phi_sub(phi_sub.length())};

    lua_pushstring(L, integrals_lualib_regkey);
    lua_gettable(L, LUA_REGISTRYINDEX);
    gsl_integration_workspace **ws =
        (gsl_integration_workspace**)lua_touserdata(L, -1);

    method_funcs homotopy_method = {pwlinear_cubic_homotopy_phi,
                                    pwlinear_cubic_homotopy_dphi,
                                    compute_convolution_s,
                                    removed_singularity_s,
                                    far_upstream_integral_s,
                                    far_downstream_integral_s};

    u_s_worker(v_s,
               phi_sub,
               dphi_sub,
               beta,
               beta_sub,
               ffus_params,
               ffds_params,
               s,
               ws,
               homotopy_method,
               u_s);
  
    veclua_pushtable(L, u_s);

    return 1;

}


/** Entry point to compute surface \f$u_s(\beta_{j+1/2})\f$ using \f$v\f$ for
    symmetric free surface problems.
  
    Arguments are passed by the lua stack
    - v_s \f$(n-1)\f$-vector; the horizontal component of velocity on the free
          surface at \f$\beta_{j+1/2}\f$.
    - phi_sub \f$(m)\f$-vector; 'clustering' locations in computational grid;
              \f$\phi^{[j]}\f$.
    - dphi_sub \f$(m)\f$-vector; value of \f$\phi'\f$ at \f$\phi^{[j]}\f$.
    - beta \f$(n)\f$-vector; the grid points \f$\beta_j\f$.
    - beta_sub \f$(m)\f$-vector; \f$\beta\f$ at \f$\phi^{[j]}\f$.
    - D scalar; coefficent of linearised far field velocity.
    - lambda scalar; exponent of linearised far field velocity.
    - s scalar; optional homotopy parameter.

 \param L Lua stack pointer.
 
*/
int u_homotopy_symmetric_s(lua_State* L) {

    static const int NARG = 8;

    const int n = veclua_veclength(L, -NARG+0);
    const int n_sub = veclua_veclength(L, -NARG+1);

    const std::unique_ptr<real[]> heap_v(new real [(2*n)+(n-1)+(3*n_sub)]);

    /* free surface variables */
    const real_vector v_s(n, &(heap_v[0]), L, -NARG+0);
    /* grid variables */
    const real_vector     beta(    n,               &(heap_v[n]), L, -NARG+3);
    const real_vector  phi_sub(n_sub,             &(heap_v[2*n]), L, -NARG+1);
    const real_vector dphi_sub(n_sub,     &(heap_v[(2*n)+n_sub]), L, -NARG+2);
    const real_vector beta_sub(n_sub, &(heap_v[(2*n)+(2*n_sub)]), L, -NARG+4);

    /* far field parameters */
    const real D      = lua_tonumber(L, -NARG+5);
    const real lambda = lua_tonumber(L, -NARG+6);

    /* homotopy parameter */
    const real s = lua_tonumber(L, -NARG+7);

    lua_pop(L, NARG);

    real_vector u_s(n, &(heap_v[(2*n)+(3*n_sub)]));

    ff_params ffus_params = {lambda, 0, D, -phi_sub(phi_sub.length())};
    ff_params ffds_params = {lambda, 0, D, phi_sub(phi_sub.length())};

    lua_pushstring(L, integrals_lualib_regkey);
    lua_gettable(L, LUA_REGISTRYINDEX);
    gsl_integration_workspace **ws =
        (gsl_integration_workspace**)lua_touserdata(L, -1);

    method_funcs homotopy_method = {pwlinear_cubic_homotopy_phi,
                                    pwlinear_cubic_homotopy_dphi,
                                    compute_convolution_symmetric_s,
                                    removed_singularity_symmetric_s,
                                    far_upstream_integral_s,
                                    far_downstream_integral_s};

    u_s_worker(v_s,
               phi_sub,
               dphi_sub,
               beta,
               beta_sub,
               ffus_params,
               ffds_params,
               s,
               ws,
               homotopy_method,
               u_s);
  
    veclua_pushtable(L, u_s);

    return 1;

}


