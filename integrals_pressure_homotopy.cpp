/** \file integrals_pressure_homotopy.cpp
  
*/

#include "integrals_pressure.h"
#include "integrals_pressure_formulation.h"

extern "C" int luaopen_integrals_pressure_homotopy(lua_State*);

/* Homotopy mapping */
void pwlinear_cubic_homotopy_phi(const vector&, const vector&, const vector&, const vector&, vector&, const real);
void pwlinear_cubic_homotopy_dphi(const vector&, const vector&, const vector&, const vector&, vector&, const real);

static int pwlinear_cubic_homotopy_phi_luaf(lua_State*);
static int pwlinear_cubic_homotopy_dphi_luaf(lua_State*);

/* General formulation */
static int x_homotopy_s(lua_State*);
static int y_homotopy_s(lua_State*);

static int u_homotopy_s(lua_State*);

/* Symmetric formulation */
static int x_homotopy_symmetric_s(lua_State*);
/* y_homotopy_symmetric_s does the same thing as y_homotopy_s */
/*static int y_homotopy_symmetric_s(lua_State*);*/

static int u_homotopy_symmetric_s(lua_State*);

extern "C" const char integrals_lualib_name[]   = "integrals_pressure";
extern "C" const char integrals_lualib_regkey[] = "integrals_pressure_rk";
extern "C" const luaL_reg integralsR[] = {
{"u_s", u_homotopy_s},
{"x_s", x_homotopy_s},
{"y_s", y_homotopy_s},
{"u_symmetric_s", u_homotopy_symmetric_s},
{"x_symmetric_s", x_homotopy_symmetric_s},
{"y_symmetric_s", y_homotopy_s}, /* used to be y_homotopy_symmetric_s */
{"phi", pwlinear_cubic_homotopy_phi_luaf},
{"dphi", pwlinear_cubic_homotopy_dphi_luaf},
{NULL, NULL}
};

int luaopen_integrals_pressure_homotopy(lua_State *L)
{
  /** \todo May want to clarify the naming of the library. */
  return luaopen_integrals_pressure(L);
}

/*
  Assumes \f$\beta\f$ is in increasing order.
*/

/** Calculates a piece-wise polynomial relationship between \f$\beta\f$ and \f$\phi\f$ based on a homotopy between a linear and a cubic interpolating spline.

  Calculates a piece-wise polynomial relationship between two variables,
  \f$\beta\f$ and \f$\phi\f$, which is based on a homotopy parameter \f$s\f$
  which varies the relationship between a piecewise linear spline and a piecewise
  cubic spline. The genereal piecewise polynomial/homotopy is defined on the unit
  intervals of \f$\beta\f$ given by \f$[i-1,i]\f$ for \f$i = 1,2,\ldots,n\f$.
  Then the piecewise polynomial/homotopy is
  \f[ p(\beta) = p_i(\beta) = (1-s) A_i(\beta) + sS_i(\beta) \f]
  for \f$\beta \in [i-1,i]\f$.

  The data required for the polynomial are $\phi_i$ and $\phi_i'$ for
  \f$i = 0, 1, \ldots, n\f$.

  \f$A_i\f$ is the piecewise linear interpolating function that satisfies two
  equations based on the data:

  \f$A_i(i) = \phi_i\f$ and \f$A_i(i-1) = \phi_{i-1}\f$ for each 
  \f$i = 1,2, \ldots, n\f$.

  \f$S_i\f$ is the piecewise cubic interpolating spline that satisfies four
  equations based on the data:

  \f$S_i(i) = \phi_i\f$, \f$S_i(i-1) = \phi_{i-1}\f$, \f$S_i'(i) = \phi_i'\f$,
  and \f$S_i'(i-1) = \phi_{i-1}'\f$ for each 
  \f$i = 1,2, \ldots, n\f$.

  When \f$s = 0\f$ clearly the function is piecewise linear, and when
  \f$s = 1\f$ it is a piecewise cubic.

  \param beta A \f$(m)\f$-vector specifying the points at which the polynomial will be evaluated. 
  \param beta_sub A \f$(n+1)\f$-vector specifying data to construct piecewise polynomial.
  \param phi_sub A \f$(n+1)\f$-vector specifying data to construct the piecewise polynomial.
  \param dphi_sub  A \f$(n+1)\f$-vector specifying data to construct the piecewise polynomial.
  \param o A \f$(m)\f$-vector with the output of the function evaluated at the input \f$\beta\f$
  \param s A homotopy parameter that should be between 0 and 1.

*/
void pwlinear_cubic_homotopy_phi(const vector&     beta,
                                 const vector& beta_sub,
                                 const vector&  phi_sub,
                                 const vector& dphi_sub,
                                 vector              &o,
                                 const real           s)
{
  unsigned int m = phi_sub.length();

  unsigned int j = 1;
  for(unsigned int i = 1; i <= beta.length(); i++) {

    if(beta(i) <= beta_sub(1)) {
      o(i) = phi_sub(1) + (beta(i) - beta_sub(1)) * dphi_sub(1);
    } else if(beta(i) > beta_sub(m)) {
      o(i) = phi_sub(m) + (beta(i) - beta_sub(m)) * dphi_sub(m);
    } else {
      while(beta_sub(j) < beta(i)) {
        j++;
      }
      
      // beta_sub(j) is now right-hand end point, beta_sub(j-1) is left-hand end point

      // Okay so the following terms help with the factors of the basis polynomials
      real dbeta = beta_sub(j) - beta_sub(j-1);

      real bi_sub_bj0 = beta(i) - beta_sub(j-1);
      real bi_sub_bj1 = beta(i) - beta_sub(j);

      real zero1 = (3.0*beta_sub(j-1) - beta_sub(j)) / 2.0;
      real zero2 = (3.0*beta_sub(j) - beta_sub(j-1)) / 2.0;

      // Calculate values of basis polynomials
      real f1 = bi_sub_bj1 * bi_sub_bj1 *
                (beta(i) - zero1) / (dbeta*dbeta*(beta_sub(j-1) - zero1));
      real f2 = bi_sub_bj0 * bi_sub_bj0 *
                (beta(i) - zero2) / (dbeta*dbeta*(beta_sub(j) - zero2));
      real f3 = bi_sub_bj0 * bi_sub_bj1 * bi_sub_bj1 / (dbeta*dbeta);
      real f4 = bi_sub_bj1 * bi_sub_bj0 * bi_sub_bj0 / (dbeta*dbeta);

      // Construct phi using homotopy between linear interpolation and the cubic
      //   interpolation with fixed, continuous, derivatives
      o(i)  = (1.0-s) * ((phi_sub(j)*bi_sub_bj0) - (phi_sub(j-1)*bi_sub_bj1) ) / (dbeta);
      o(i) += s*(f1*phi_sub(j-1)  + f2*phi_sub(j) +
                 f3*dphi_sub(j-1) + f4*dphi_sub(j));
    }
  }
}

void pwlinear_cubic_homotopy_dphi(const vector&     beta,
                                  const vector& beta_sub,
                                  const vector&  phi_sub,
                                  const vector& dphi_sub,
                                  vector              &o,
                                  const real           s)
{
  unsigned int m = phi_sub.length();

  unsigned int j = 1;
  for(unsigned int i = 1; i <= beta.length(); i++) {

    if(beta(i) <= beta_sub(1)) {
      o(i) = dphi_sub(1);
    } else if(beta(i) > beta_sub(m)) {
      o(i) = dphi_sub(m);
    } else {
      while(beta_sub(j) < beta(i)) {
        j++;
      }

      // beta_sub(j) is now right-hand end point, beta_sub(j-1) is left-hand end point

      // Okay so the following terms help with the factors of the basis polynomials
      real dbeta = beta_sub(j) - beta_sub(j-1);

      real bi_sub_bj0 = beta(i) - beta_sub(j-1);
      real bi_sub_bj1 = beta(i) - beta_sub(j);

      real zero1 = (3.0*beta_sub(j-1) - beta_sub(j)) / 2.0;
      real zero2 = (3.0*beta_sub(j) - beta_sub(j-1)) / 2.0;

      // Calculate values of basis polynomials
      real df1 = bi_sub_bj1 * (2.0*(beta(i) - zero1) + bi_sub_bj1) /
                 (dbeta*dbeta*(beta_sub(j-1) - zero1));
      real df2 = bi_sub_bj0 * (2.0*(beta(i) - zero2) + bi_sub_bj0) /
                 (dbeta*dbeta*(beta_sub(j) - zero2));
      real df3 = bi_sub_bj1*(2.0*bi_sub_bj0 + bi_sub_bj1) / (dbeta*dbeta);
      real df4 = bi_sub_bj0*(2.0*bi_sub_bj1 + bi_sub_bj0) / (dbeta*dbeta);

      // Construct phi using homotopy between linear interpolation and the cubic
      //   interpolation with fixed, continuous, derivatives
      o(i)  = (1.0-s) * (phi_sub(j) - phi_sub(j-1)) / (dbeta);
      o(i) += s*(df1*phi_sub(j-1)  + df2*phi_sub(j) +
                 df3*dphi_sub(j-1) + df4*dphi_sub(j));
    }
  }
}

int pwlinear_cubic_homotopy_phi_luaf(lua_State* L)
{
  static const int NARG = 5;

  vector beta     = veclua_tovector(L, -NARG+0);
  vector beta_sub = veclua_tovector(L, -NARG+1);
  vector phi_sub  = veclua_tovector(L, -NARG+2);
  vector dphi_sub = veclua_tovector(L, -NARG+3);
  real   s        = lua_tonumber(L, -NARG+4);

  lua_pop(L, NARG);

  unsigned int m = phi_sub.length();
  unsigned int j = 1;

  lua_newtable(L);

  for(unsigned int i = 1; i <= beta.length(); i++) {
    real o_i = 0.0;

    if(beta(i) <= beta_sub(1)) {
      o_i = phi_sub(1) + (beta(i) - beta_sub(1)) * dphi_sub(1);
    } else if(beta(i) > beta_sub(m)) {
      o_i = phi_sub(m) + (beta(i) - beta_sub(m)) * dphi_sub(m);
    } else {
      while(beta_sub(j) < beta(i)) {
        j++;
      }
      // beta_sub(j) is now right-hand end point, beta_sub(j-1) is left-hand end point

      // Okay so the following terms help with the factors of the basis polynomials
      real dbeta = beta_sub(j) - beta_sub(j-1);

      real bi_sub_bj0 = beta(i) - beta_sub(j-1);
      real bi_sub_bj1 = beta(i) - beta_sub(j);

      real zero1 = (3.0*beta_sub(j-1) - beta_sub(j)) / 2.0;
      real zero2 = (3.0*beta_sub(j) - beta_sub(j-1)) / 2.0;

      // Calculate values of basis polynomials
      real f1 = bi_sub_bj1 * bi_sub_bj1 *
                (beta(i) - zero1) / (dbeta*dbeta*(beta_sub(j-1) - zero1));
      real f2 = bi_sub_bj0 * bi_sub_bj0 *
                (beta(i) - zero2) / (dbeta*dbeta*(beta_sub(j) - zero2));
      real f3 = bi_sub_bj0 * bi_sub_bj1 * bi_sub_bj1 / (dbeta*dbeta);
      real f4 = bi_sub_bj1 * bi_sub_bj0 * bi_sub_bj0 / (dbeta*dbeta);

      // Construct phi using homotopy between linear interpolation and the cubic
      //   interpolation with fixed, continuous, derivatives
      o_i  = (1.0-s) * ((phi_sub(j)*bi_sub_bj0) - (phi_sub(j-1)*bi_sub_bj1) ) / (dbeta);
      o_i += s*(f1*phi_sub(j-1)  + f2*phi_sub(j) +
                f3*dphi_sub(j-1) + f4*dphi_sub(j));
    }

    lua_pushnumber(L, i);
    lua_pushnumber(L, o_i);
    lua_rawset(L, -3);
  }

  return 1;
}


int pwlinear_cubic_homotopy_dphi_luaf(lua_State* L)
{
  static const int NARG = 5;

  vector beta     = veclua_tovector(L, -NARG+0);
  vector beta_sub = veclua_tovector(L, -NARG+1);
  vector phi_sub  = veclua_tovector(L, -NARG+2);
  vector dphi_sub = veclua_tovector(L, -NARG+3);
  real   s        = lua_tonumber(L, -NARG+4);

  lua_pop(L, NARG);

  unsigned int m = phi_sub.length();
  unsigned int j = 1;

  lua_newtable(L);

  for(unsigned int i = 1; i <= beta.length(); i++) {
    real o_i = 0.0;

    if(beta(i) <= beta_sub(1)) {
      o_i = dphi_sub(1);
    } else if(beta(i) > beta_sub(m)) {
      o_i = dphi_sub(m);
    } else {
      while(beta_sub(j) < beta(i)) {
        j++;
      }

      // beta_sub(j) is now right-hand end point, beta_sub(j-1) is left-hand end point

      // Okay so the following terms help with the factors of the basis polynomials
      real dbeta = beta_sub(j) - beta_sub(j-1);

      real bi_sub_bj0 = beta(i) - beta_sub(j-1);
      real bi_sub_bj1 = beta(i) - beta_sub(j);

      real zero1 = (3.0*beta_sub(j-1) - beta_sub(j)) / 2.0;
      real zero2 = (3.0*beta_sub(j) - beta_sub(j-1)) / 2.0;

      // Calculate values of basis polynomials
      real df1 = bi_sub_bj1 * (2.0*(beta(i) - zero1) + bi_sub_bj1) /
                 (dbeta*dbeta*(beta_sub(j-1) - zero1));
      real df2 = bi_sub_bj0 * (2.0*(beta(i) - zero2) + bi_sub_bj0) /
                 (dbeta*dbeta*(beta_sub(j) - zero2));
      real df3 = bi_sub_bj1*(2.0*bi_sub_bj0 + bi_sub_bj1) / (dbeta*dbeta);
      real df4 = bi_sub_bj0*(2.0*bi_sub_bj1 + bi_sub_bj0) / (dbeta*dbeta);

      // Construct phi using homotopy between linear interpolation and the cubic
      //   interpolation with fixed, continuous, derivatives
      o_i  = (1.0-s) * (phi_sub(j) - phi_sub(j-1)) / (dbeta);
      o_i += s*(df1*phi_sub(j-1)  + df2*phi_sub(j) +
                df3*dphi_sub(j-1) + df4*dphi_sub(j));
    }

    lua_pushnumber(L, i);
    lua_pushnumber(L, o_i);
    lua_rawset(L, -3);
  }

  return 1;
}


/** Entry point to compute surface \f$x_s(\beta_{j+1/2})\f$ using \f$u\f$ and \f$v\f$.

   Arguments are passed by the Lua stack
 - u_mid \f$(n-1)\f$-vector; the horizontal component of velocity on the free surface at \f$\beta_{j+1/2}\f$.
 - v_mid \f$(n-1)\f$-vector; the vertical component of velocity on the free surface at \f$\beta_{j+1/2}\f$.
 - beta \f$(n)\f$-vector; the grid points \f$\beta_j\f$.
 - beta_sub \f$(m)\f$-vector; \f$\beta\f$ at \f$\phi^{[j]}\f$.
 - phi_sub \f$(m)\f$-vector; 'clustering' locations in computational grid; \f$\phi^{[j]}\f$.
 - dphi_sub \f$(m)\f$-vector; value of \f$\phi'\f$ at \f$\phi^{[j]}\f$.
 - s scalar; optional homotopy parameter.
 \param L Lua stack pointer.

*/
int x_homotopy_s(lua_State* L)
{
  static const int NARG = 7;
  
  vector u_mid    = veclua_tovector(L, -NARG+0);
  vector v_mid    = veclua_tovector(L, -NARG+1);
  vector beta     = veclua_tovector(L, -NARG+2);
  vector beta_sub = veclua_tovector(L, -NARG+3);
  vector phi_sub  = veclua_tovector(L, -NARG+4);
  vector dphi_sub = veclua_tovector(L, -NARG+5);
  real   s        = lua_tonumber(L, -NARG+6);

  real x0;

  lua_pop(L, NARG);
  
  vector x_s = x_s_worker(u_mid,
                          v_mid,
                           beta,
                       beta_sub,
                        phi_sub,
                       dphi_sub,
                              s,
   pwlinear_cubic_homotopy_dphi,
                             x0);

  lua_newtable(L);
  for(unsigned int j = 1; j <= x_s.length(); j++) {
    lua_pushnumber(L, j);
    lua_pushnumber(L, x_s(j));
    lua_rawset(L, -3);
  }

  lua_pushnumber(L, x0);

  return 2;
}

/** Entry point to compute surface \f$y_s(\beta_{j+1/2})\f$ using \f$u\f$ and \f$v\f$.
  
   Arguments are passed by the lua stack
 - u_mid \f$(n-1)\f$-vector; the horizontal component of velocity on the free surface at \f$\beta_{j+1/2}\f$.
 - v_mid \f$(n-1)\f$-vector; the vertical component of velocity on the free surface at \f$\beta_{j+1/2}\f$.
 - beta \f$(n)\f$-vector; the grid points \f$\beta_j\f$.
 - beta_sub \f$(m)\f$-vector; \f$\beta\f$ at \f$\phi^{[j]}\f$.
 - phi_sub \f$(m)\f$-vector; 'clustering' locations in computational grid; \f$\phi^{[j]}\f$.
 - dphi_sub \f$(m)\f$-vector; value of \f$\phi'\f$ at \f$\phi^{[j]}\f$.
 - s scalar; optional homotopy parameter.
 \param L Lua stack pointer.
*/

int y_homotopy_s(lua_State* L)
{
  static const int NARG = 7;
  
  vector u_mid     = veclua_tovector(L, -NARG+0);
  vector v_mid     = veclua_tovector(L, -NARG+1);
  vector beta      = veclua_tovector(L, -NARG+2);
  vector beta_sub  = veclua_tovector(L, -NARG+3);
  vector phi_sub   = veclua_tovector(L, -NARG+4);
  vector dphi_sub  = veclua_tovector(L, -NARG+5);
  real   s         = lua_tonumber(L, -NARG+6);

  real eta0;
 
  lua_pop(L, NARG);

  vector y_s = y_s_worker(u_mid, 
                          v_mid,
                           beta, 
                       beta_sub,
                        phi_sub,
                       dphi_sub,
                              s,
   pwlinear_cubic_homotopy_dphi,
                           eta0);

  lua_newtable(L);
  for(unsigned int j = 1; j <= y_s.length(); j++) {
    lua_pushnumber(L, j);
    lua_pushnumber(L, y_s(j));
    lua_rawset(L, -3);
  }

  lua_pushnumber(L, eta0);
  
  return 2;
}

/* Entry point to compute surface \f$u_s(\beta_{j+1/2})\f$ using \f$v\f$.
  
   Arguments are passed by the lua stack
 - v_s \f$(n-1)\f$-vector; the horizontal component of velocity on the free surface at \f$\beta_{j+1/2}\f$.
 - phi_sub \f$(m)\f$-vector; 'clustering' locations in computational grid; \f$\phi^{[j]}\f$.
 - dphi_sub \f$(m)\f$-vector; value of \f$\phi'\f$ at \f$\phi^{[j]}\f$.
 - beta \f$(n)\f$-vector; the grid points \f$\beta_j\f$.
 - beta_sub \f$(m)\f$-vector; \f$\beta\f$ at \f$\phi^{[j]}\f$.
 - Dv_dstream scalar; coefficent of linearised downstream velocity solution.
 - Dv_ustream scalar; coefficent of linearised upstream velocity solution.
 - lambda scalar; exponent of linearised farstream velocity solution.
 - s scalar; optional homotopy parameter.
 \param L Lua stack pointer.  */
int u_homotopy_s(lua_State* L)
{
  static const int NARG = 9;

  const vector v_s      = veclua_tovector(L, -NARG+0);
  const vector phi_sub  = veclua_tovector(L, -NARG+1);
  const vector dphi_sub = veclua_tovector(L, -NARG+2);
  const vector beta     = veclua_tovector(L, -NARG+3);
  const vector beta_sub = veclua_tovector(L, -NARG+4);
  const real Dv_dstream = lua_tonumber(L, -NARG+5);
  const real Dv_ustream = lua_tonumber(L, -NARG+6);
  const real lambda     = lua_tonumber(L, -NARG+7);
  const real s          = lua_tonumber(L, -NARG+8);

  lua_pop(L, NARG);

  linear_params lin_us = {lambda, Dv_ustream, phi_sub(1)};
  linear_params lin_ds = {lambda, Dv_dstream, phi_sub(phi_sub.length())};

  lua_pushstring(L, integrals_lualib_regkey);
  lua_gettable(L, LUA_REGISTRYINDEX);
  gsl_integration_workspace **ws = (gsl_integration_workspace**)lua_touserdata(L, -1);

  method_funcs homotopy_method = {&pwlinear_cubic_homotopy_phi,
                                 &pwlinear_cubic_homotopy_dphi,
                                        &compute_convolution_s,
                                &removed_singularity_general_s,
                              &compute_far_upstream_integral_s,
                            &compute_far_downstream_integral_s};

  const vector u_s = u_s_worker(v_s,
                            phi_sub,
                           dphi_sub,
                               beta,
                           beta_sub,
                             lin_us,
                             lin_ds,
                                  s,
                                 ws,
                    homotopy_method);
  
    
  lua_newtable(L);

  for(unsigned int i = 1; i <= u_s.length(); i++) {
    lua_pushnumber(L, i);
    lua_pushnumber(L, u_s(i));
    lua_rawset(L, -3);
  }

  return 1;
}

/** Entry point to compute surface \f$x_s(\beta_{j+1/2})\f$ using \f$u\f$ and \f$v\f$.
  
   Arguments are passed by the lua stack
 - u_mid \f$(n-1)\f$-vector; the horizontal component of velocity on the free surface at \f$\beta_{j+1/2}\f$.
 - v_mid \f$(n-1)\f$-vector; the vertical component of velocity on the free surface at \f$\beta_{j+1/2}\f$.
 - beta \f$(n)\f$-vector; the grid points \f$\beta_j\f$.
 - beta_sub \f$(m)\f$-vector; \f$\beta\f$ at \f$\phi^{[j]}\f$.
 - phi_sub \f$(m)\f$-vector; 'clustering' locations in computational grid; \f$\phi^{[j]}\f$.
 - dphi_sub \f$(m)\f$-vector; value of \f$\phi'\f$ at \f$\phi^{[j]}\f$.
 - s scalar; optional homotopy parameter.
 \param L Lua stack pointer.
*/
int x_homotopy_symmetric_s(lua_State* L)
{
  static const int NARG = 7;
  
  vector u_mid    = veclua_tovector(L, -NARG+0);
  vector v_mid    = veclua_tovector(L, -NARG+1);
  vector beta     = veclua_tovector(L, -NARG+2);
  vector beta_sub = veclua_tovector(L, -NARG+3);
  vector phi_sub  = veclua_tovector(L, -NARG+4);
  vector dphi_sub = veclua_tovector(L, -NARG+5);
  real   s        = lua_tonumber(L, -NARG+6);

  real x0;

  lua_pop(L, NARG);
  
  vector x_s = x_s_worker(u_mid,
                          v_mid,
                           beta,
                       beta_sub,
                        phi_sub,
                       dphi_sub,
                              s,
   pwlinear_cubic_homotopy_dphi,
                             x0);

  lua_newtable(L);
  for(unsigned int j = 1; j <= x_s.length(); j++) {
    lua_pushnumber(L, j);
    lua_pushnumber(L, x_s(j)-x0);
    lua_rawset(L, -3);
  }

  return 1;
}


/* \todo Template for u_homotopy_s ... needs some thought still */
int u_homotopy_symmetric_s(lua_State* L)
{
  static const int NARG = 8;

  const vector v_s      = veclua_tovector(L, -NARG+0);
  const vector phi_sub  = veclua_tovector(L, -NARG+1);
  const vector dphi_sub = veclua_tovector(L, -NARG+2);
  const vector beta     = veclua_tovector(L, -NARG+3);
  const vector beta_sub = veclua_tovector(L, -NARG+4);
  const real Dv         = lua_tonumber(L, -NARG+5);
  const real lambda     = lua_tonumber(L, -NARG+6);
  const real s          = lua_tonumber(L, -NARG+7);

  lua_pop(L, NARG);

  linear_params lin_us = {lambda, Dv, -phi_sub(phi_sub.length())};
  linear_params lin_ds = {lambda, Dv,  phi_sub(phi_sub.length())};

  lua_pushstring(L, integrals_lualib_regkey);
  lua_gettable(L, LUA_REGISTRYINDEX);
  gsl_integration_workspace **ws = (gsl_integration_workspace**)lua_touserdata(L, -1);

  method_funcs homotopy_method = {&pwlinear_cubic_homotopy_phi,
                                 &pwlinear_cubic_homotopy_dphi,
                              &compute_convolution_symmetric_s,
                              &removed_singularity_symmetric_s,
                              &compute_far_upstream_integral_s,
                            &compute_far_downstream_integral_s};

  const vector u_s = u_s_worker(v_s,
                            phi_sub,
                           dphi_sub,
                               beta,
                           beta_sub,
                             lin_us,
                             lin_ds,
                                  s,
                                 ws,
                    homotopy_method);
  
    
  lua_newtable(L);

  for(unsigned int i = 1; i <= u_s.length(); i++) {
    lua_pushnumber(L, i);
    lua_pushnumber(L, u_s(i));
    lua_rawset(L, -3);
  }

  return 1;
}
