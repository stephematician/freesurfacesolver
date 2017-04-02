/** \file integrals_pressure.cpp

  C implementation of integral formulas.

*/

#include "integrals_pressure.h"

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
 
  /** \todo Fix which name to use for the gsl_integration_workspace registry 'key' */
  lua_pushstring(L, integrals_lualib_regkey);
  // Make new userdata at the top of the stack
  gsl_integration_workspace **ws = (gsl_integration_workspace**)lua_newuserdata(L, sizeof(gsl_integration_workspace **));
  *ws = gsl_integration_workspace_alloc(1000);

  // Make a new metatable
  lua_newtable(L);
  // Garbage collection method
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, integrals_pressure_clean_up);
  lua_settable(L, -3);
  // Set metatable
  lua_setmetatable(L, -2);
  
  // Should now have stack (T->B) -> userdata, "integrals_box", luaL register table
  lua_settable(L, LUA_REGISTRYINDEX);
  
 // Should now have stack (T->B) -> luaL register table  
  return 1;
}

/* Clean up function for the gsl_integration_workspace */
int integrals_pressure_clean_up(lua_State* L)
{
  lua_pushstring(L, integrals_lualib_regkey);
  lua_gettable(L, LUA_REGISTRYINDEX);
  gsl_integration_workspace **ws = (gsl_integration_workspace**)lua_touserdata(L, -1);
  
  gsl_integration_workspace_free(*ws);

  return 0;
}


/** Computes \f$u_s(\beta)\f$ at mid points, given flat topography.

  \param v_s A \f$(n+1)\f$-vector containing the vertical component of fluid velocity at grid points. 
  \param phi_sub A \f$(n+1)\f$-vector specifying data to construct the piecewise polynomial.
  \param dphi_sub  A \f$(n+1)\f$-vector specifying data to construct the piecewise polynomial.
  \param beta A \f$(n)\f$-vector specifying the location in \f$\beta\f$ of the desired mid-points.
  \param beta_sub A \f$(M+1)\f$-vector specifying data to construct piecewise polynomial.
  \param lin_ustream Parameters for calculating the upstream linearised integral.
  \param lin_dstream Parameters for calculating the downstream linearised integral.
  \param s A homotopy parameter that should be between 0 and 1.
  \param ws Workspace for GSL integration used/instantiated by this package.
  \param method Structure containing functions which help define the problem (e.g. symmetric)
 
  \todo Comment on what this actually calculates (check note5, note6 and thesis?).
*/


vector u_s_worker(const vector&                v_s,
                  const vector&            phi_sub,
                  const vector&           dphi_sub,
                  const vector&               beta,
                  const vector&           beta_sub,
                  const linear_params& lin_ustream,
                  const linear_params& lin_dstream,
                  const real                     s,
                  gsl_integration_workspace**   ws,
                  method_funcs&             method)
{

  unsigned int n = beta.length();

  vector phi  = vector(n);
  vector dphi = vector(n);

  vector integrand = vector(n);

  vector beta_mid = vector(n-1);
  vector phi_mid  = vector(n-1);

  vector u_s = vector(n-1);

  for(unsigned int i = 1; i <= n-1; i++) {
    beta_mid(i) = (beta(i) + beta(i+1)) / 2.0;
  }

  (*(method.transform_phi))(beta, beta_sub, phi_sub, dphi_sub, phi, s);
  (*(method.transform_dphi))(beta, beta_sub, phi_sub, dphi_sub, dphi, s);
  
  (*(method.transform_phi))(beta_mid, beta_sub, phi_sub, dphi_sub, phi_mid, s);
  
  vector removed_singularity = vector(n-1);
  (*(method.removed_singularity))(v_s, phi, phi_mid, removed_singularity);

  vector downstream_term = vector(n-1);
  (*(method.linear_downstream_int))(phi_mid,
                                    lin_dstream.phi_match,
                                    lin_dstream.D,
                                    lin_dstream.lambda,
                                    downstream_term);

  vector upstream_term = vector(n-1);
  (*(method.linear_upstream_int))(phi_mid,
                                  lin_ustream.phi_match,
                                  lin_ustream.D,
                                  lin_ustream.lambda,
                                  upstream_term);

  for(unsigned int i = 1; i <= n-1; i++) {
    /* From the formulation, begin with u = 1.0 */
    u_s(i) = 1.0;

    (*(method.convolution_terms))(i, v_s, phi, dphi, phi_mid, integrand);

    /* This implementation of a 'trapezoidal' type integral has been tested */
    for(unsigned int j = 1; j <= n; j++) {
      real delta = (j == 1 ? beta(2) - beta(1) :
                             (j == n ? beta(n) - beta(n-1) :
                                       beta(j+1) - beta(j-1)));
      u_s(i) += delta * integrand(j) / 2.0;
    }

    /* Expression for removal of singularity */
    u_s(i) += removed_singularity(i);

    /* Contribution to integral from far downstream */
    u_s(i) += downstream_term(i);

    /* Contribution to integral from far upstream */
    u_s(i) += upstream_term(i);
  }

  return u_s;
}

vector y_s_worker(const vector           &u_mid,
                  const vector           &v_mid,
                  const vector            &beta,
                  const vector        &beta_sub,
                  const vector         &phi_sub,
                  const vector        &dphi_sub,
                  const real                  s,
                  transform_func transform_dphi,
                  real                    &eta0)
{
  unsigned int n = beta.length();
    
  vector beta_mid = vector(n-1);
  for(unsigned int i = 1; i <= n-1; i++) {
    beta_mid(i) = (beta(i) + beta(i+1)) / 2.0;
  }
  
  vector dphi_mid = vector(n-1);
  (*transform_dphi)(beta_mid, beta_sub, phi_sub, dphi_sub, dphi_mid, s);
  
  vector y_s = vector(n-1);
  
  real f_j = -dphi_mid(n-1) * v_mid(n-1) / \
             (u_mid(n-1)*u_mid(n-1) + v_mid(n-1)*v_mid(n-1));
  real M_j = f_j * (beta(n) - beta(n-1));
  y_s(n-1) = M_j / 2.0;

  for(unsigned int j = n-2; j >= 1; j--) {
    f_j = -dphi_mid(j) * v_mid(j) / (u_mid(j)*u_mid(j) + v_mid(j)*v_mid(j));
    M_j = M_j + f_j * (beta(j+1) - beta(j));
    
    y_s(j) = M_j - f_j * (beta(j+1) - beta(j)) / 2.0;
  }
  
  eta0 = M_j;
  
  return y_s;
}


vector x_s_worker(const vector           &u_mid,
                  const vector           &v_mid,
                  const vector            &beta,
                  const vector        &beta_sub,
                  const vector         &phi_sub,
                  const vector        &dphi_sub,
                  const real                  s,
                  transform_func transform_dphi,
                  real                      &x0)
{
  
  unsigned int n = beta.length();
    
  vector beta_mid = vector(n-1);
  for(unsigned int i = 1; i <= n-1; i++) {
    beta_mid(i) = (beta(i) + beta(i+1)) / 2.0;
  }
  
  vector dphi_mid = vector(n-1);
  (*transform_dphi)(beta_mid, beta_sub, phi_sub, dphi_sub, dphi_mid, s);
  
  vector x_s = vector(n-1);
  
  real f_j = -dphi_mid(n-1) * u_mid(n-1) / \
             (u_mid(n-1)*u_mid(n-1) + v_mid(n-1)*v_mid(n-1)); 
  real M_j = f_j * (beta(n) - beta(n-1));
  x_s(n-1) = M_j / 2.0;

  for(unsigned int j = n-2; j >= 1; j--) {
    f_j = -dphi_mid(j) * u_mid(j) / (u_mid(j)*u_mid(j) + v_mid(j)*v_mid(j));
    M_j = M_j + f_j * (beta(j+1) - beta(j));
    
    x_s(j) = M_j - f_j * (beta(j+1) - beta(j)) / 2.0;
  }

  x0 = M_j;

  return x_s;
}
