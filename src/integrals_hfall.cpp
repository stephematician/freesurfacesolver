/** \file integrals_hfall.cpp

  C implementation of integral formulas.

*/

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "vector.h"

#define pi M_PI

extern "C" int luaopen_integrals_hfall(lua_State *);
static int tau_f(lua_State*);
static int tau_s(lua_State*);
static int tau_b(lua_State*);

static int y_f(lua_State*);
static int y_b(lua_State*);
static int x_f(lua_State*);
static int x_b(lua_State*);

static const luaL_reg integralsR[] = {
{"tau", tau_f},
{"tau_s", tau_s},
{"tau_b", tau_b},
{"y", y_f},
{"y_b", y_b},
{"x", x_f},
{"x_b", x_b},
{NULL, NULL}
};

int luaopen_integrals_hfall(lua_State *L) {
  luaL_register(L, "integrals", integralsR);
  return 1;
}

/** Computes \f$\tau(\phi)\f$ at mid points, given unperturbed topography

 Computes the logarithm of the speed on the free surface given by
 \f$\tau(\phi_{i+1/2}) = \tau((\phi_i + \phi_{i+1})/2)\f$ using the formula
 \f[ \tau(\phi_{i+1/2}) = \sum_{j=1}^{N} w_j \frac{\theta(\phi_j) 
                      e^{\pi\phi_{j}}}{e^{\pi\phi_{j}} - e^{\pi\phi_{j+1/2}}}
                      \f]
                      
 Arguments are passed by the lua stack
 - theta, a table of \f$N\f$ values of the angle of the free surface,
 - phi, a table of \f$N\f$ value of phi (the computational grid),
 - exp_phi, a table representing \f$e^{\pi\phi}\f$,
 - delta, a real number, the grid spacing.
 
 
 \param L Lua stack pointer
*/
int tau_f(lua_State* L)
{
  static const int NARG = 4;
  
  vector theta   = veclua_tovector(L, -NARG+0);
  vector phi     = veclua_tovector(L, -NARG+1);
  vector exp_phi = veclua_tovector(L, -NARG+2);
  real delta     = lua_tonumber(L, -NARG+3);
  lua_pop(L, NARG);
  
  unsigned int n = theta.length();
  lua_newtable(L);

  for(unsigned int i = 1; i <= n-1; i++) {
    real tau_i = 0;
    /* Note : phi_ihalf = sqrt(exp_phi(i) * exp_phi(i+1)) produces bad results
       even though it is mathematically equivalent */
    real phi_ihalf = exp(pi * (phi(i) + phi(i+1)) / 2.0);
    
    for(unsigned int j = 1; j <= n; j++) {
      real w = (j == 1 || j == n) ? (delta / 2.0) : delta;
      tau_i += w * theta(j) * exp_phi(j) / (exp_phi(j) - phi_ihalf);
    }
    lua_pushnumber(L, i);
    lua_pushnumber(L, tau_i);
    lua_rawset(L, -3);
  }
  return 1;
}

/** Computes \f$\tau_s(\phi)\f$ at mid points, given triangular topography

 Computes the logarithm of the speed on the free surface given by
 \f$\tau_s(\phi_{i+1/2}) = \tau_s((\phi_i + \phi_{i+1})/2)\f$ using the formula
 \f[ \tau_s(\phi_{i+1/2}) = - \frac{\sigma}{\pi} \ln \frac{(1 + e^{\pi\phi_{i+1/2}})^2}{ 
                           |e^{\pi\phi_{\mathrm{up}}} + e^{\pi\phi_{i+1/2}}|
                           |e^{\pi\phi_{\mathrm{do}}} + e^{\pi\phi_{i+1/2}}|} +
                      \sum_{j=1}^{N} w_j \frac{\theta(\phi_j) 
                      e^{\pi\phi_{j}}}{e^{\pi\phi_{j}} - e^{\pi\phi_{j+1/2}}}
                      \f]
                      
 Arguments are passed by the lua stack
 - theta, a table of \f$N\f$ values of the angle of the free surface,
 - phi, a table of \f$N\f$ value of phi (the computational grid),
 - exp_phi, a table representing \f$e^{\pi\phi}\f$,
 - phi_up, a real number \f$\phi_{\mathrm{up}}\f$ for triangle geometry,
 - phi_down, a real number \f$\phi_{\mathrm{do}}\f$ for triangle geometry,
 - sigma, a real number \f$\sigma\f$, interior angle at base of triangle,
 - delta, a real number, the grid spacing.
 
 \param L Lua stack pointer
*/
int tau_s(lua_State* L)
{
  static const int NARG = 7;
  
  vector theta   = veclua_tovector(L, -NARG+0);
  vector phi     = veclua_tovector(L, -NARG+1);
  vector exp_phi = veclua_tovector(L, -NARG+2);
  
  real phi_up    = exp(pi * lua_tonumber(L, -NARG+3));
  real phi_down  = exp(pi * lua_tonumber(L, -NARG+4));
  real sigma     = lua_tonumber(L, -NARG+5);
  real delta     = lua_tonumber(L, -NARG+6);
  lua_pop(L, NARG);
  
  unsigned int n = theta.length();
  lua_newtable(L);
  
  for(unsigned int i = 1; i <= n-1; i++) {
    real tau_i = 0;
    /* Note : phi_ihalf = sqrt(exp_phi(i) * exp_phi(i+1)) produces bad results
       even though it is mathematically equivalent */
    real phi_ihalf = exp(pi * (phi(i) + phi(i+1)) / 2.0);
    
    for(unsigned int j = 1; j <= n; j++) {
      real w = (j == 1 || j == n) ? (delta / 2.0) : delta;
      tau_i += w * theta(j) * exp_phi(j) / (exp_phi(j) - phi_ihalf);

    }
    tau_i -= sigma * log((1 + phi_ihalf) * (1 + phi_ihalf) /
                        fabs((phi_up + phi_ihalf)*(phi_down + phi_ihalf))) / pi;
    lua_pushnumber(L, i);
    lua_pushnumber(L, tau_i);
    lua_rawset(L, -3);
  }
  return 1;
}

/** Computes \f$\tau_b(\phi)\f$ at mid points, given triangular topography

 Computes the logarithm of the speed on the channel bottom given by
 \f$\tau_b(\phi_{i+1/2}) = \tau_b((\phi_i + \phi_{i+1})/2)\f$ using the formula
 \f[ \tau_b(\phi_b) = - \frac{\sigma}{\pi} \ln \frac{(1 - e^{\pi\phi_b})^2}{ 
                           |e^{\pi\phi_{\mathrm{up}}} - e^{\pi\phi_b}|
                           |e^{\pi\phi_{\mathrm{do}}} - e^{\pi\phi_b}|} +
                      \sum_{j=1}^{N} w_j \frac{\theta(\phi_j) 
                      e^{\pi\phi_{j}}}{e^{\pi\phi_{j}} + e^{\pi\phi_b}}
                      \f]
                      
 Here \f$\phi_b\f$ is where the velocity will be computed
                      
 Arguments are passed by the lua stack
 - theta, a table of \f$N\f$ values of the angle of the free surface,
 - phi, a table of \f$N\f$ value of phi (the computational grid),
 - exp_phi, a table representing \f$e^{\pi\phi}\f$,
 - phi_b, a table of where to solve for,
 - phi_up, a real number \f$\phi_{\mathrm{up}}\f$ for triangle geometry,
 - phi_down, a real number \f$\phi_{\mathrm{do}}\f$ for triangle geometry,
 - sigma, a real number \f$\sigma\f$, interior angle at base of triangle,
 - delta, a real number, the grid spacing.
 
 \param L Lua stack pointer
*/
int tau_b(lua_State* L)
{
  static const int NARG = 8;
  
  vector theta   = veclua_tovector(L, -NARG+0);
  vector phi     = veclua_tovector(L, -NARG+1);
  vector exp_phi = veclua_tovector(L, -NARG+2);
  vector phi_b   = veclua_tovector(L, -NARG+3);
  
  real phi_up    = exp(pi * lua_tonumber(L, -NARG+4));
  real phi_down  = exp(pi * lua_tonumber(L, -NARG+5));
  real sigma     = lua_tonumber(L, -NARG+6);
  real delta     = lua_tonumber(L, -NARG+7);
  lua_pop(L, NARG);
  
  unsigned int n = theta.length();
  unsigned int m = phi_b.length();
  lua_newtable(L);

  /* working on phi_b stuff right now */
  
  for(unsigned int i = 1; i <= m-1; i++) {
    real tau_i = 0;
    /* Note : phi_ihalf = sqrt(exp_phi(i) * exp_phi(i+1)) produces bad results
       even though it is mathematically equivalent */
    real phi_b_ihalf = exp(pi * (phi_b(i) + phi_b(i+1)) / 2.0);
    
    for(unsigned int j = 1; j <= n; j++) {
      real w = (j == 1 || j == n) ? (delta / 2.0) : delta;
      tau_i += w * theta(j) * exp_phi(j) / (exp_phi(j) + phi_b_ihalf);

    }

    tau_i -= sigma * log((1 - phi_b_ihalf) * (1 - phi_b_ihalf) /
                        fabs((phi_up - phi_b_ihalf)*(phi_down - phi_b_ihalf))) / pi;
    lua_pushnumber(L, i);
    lua_pushnumber(L, tau_i);
    lua_rawset(L, -3);
  }
  return 1;
}

/** Computes free surface \f$y(\phi)\f$ at grid points, given unperturbed topography
  
  Uses a mid-point integration rule to evaluate \f$y(\phi_i)\f$,
  \f[  y(\phi_i) = y(\phi_{i+1}) - w e^{\tau_{\phi_{i+1/2}}} \sin(\theta_{i+1/2})  \f]
  with
  \f[ y(\phi_N) = 1 \f]
  
   Arguments are passed by the lua stack
 - theta, a table of \f$N\f$ values of the angle of the free surface,
 - tau, a table of \f$N-1\f$ values of tau (the computational grid),
 - delta, a real number, the grid spacing.
 
 \param L Lua stack pointer
*/
int y_f(lua_State* L)
{
  static const int NARG = 3;
  vector theta = veclua_tovector(L, -NARG+0);
  vector tau   = veclua_tovector(L, -NARG+1);
  real   delta = lua_tonumber(L, -NARG+2);
  
  lua_pop(L, NARG); 
  
  unsigned int n = theta.length();
  real y_ip1 = 1;
  
  lua_newtable(L);
  
  lua_pushnumber(L, n);
  lua_pushnumber(L, y_ip1);
  lua_rawset(L, -3);

  for(unsigned int i = n-1; i >= 1; i--) {

    real y_i = y_ip1 - delta * (exp(-tau(i)) *
                                sin((theta(i) + theta(i+1)) / 2.0));
    lua_pushnumber(L, i);
    lua_pushnumber(L, y_i);
    lua_rawset(L, -3);
    
    y_ip1 = y_i;
  }
  
  return 1;
}

/** Computes the bottom of the channel
  
  Uses a mid-point integration rule to evaluate \f$y(\phi)\f$,
  \f[  y(\phi) = \sin(\sigma) \sum_{i=1}^{N} w_i e^{-\tau_b(\phi_{i+1/2})} \f]
  with
  
   Arguments are passed by the lua stack
 - tau_b, a table of \f$N-1\f$ values of tau at mid-points,
 - sigma, a real number \f$\sigma\f$, interior angle at base of triangle,
 - delta, a real number, the grid spacing.
 
 \param L Lua stack pointer
*/
int y_b(lua_State* L)
{
  static const int NARG = 3;
  vector tau_b = veclua_tovector(L, -NARG+0);
  vector phi_b = veclua_tovector(L, -NARG+1);
  
  real sigma = lua_tonumber(L, -NARG+2);
  
  lua_pop(L, NARG); 
  
  unsigned int n = phi_b.length();
  real y_result = 0.0;
  
  for(unsigned int i = 1; i <= (n-1); i++) {
    y_result += sin(sigma) * (phi_b(i+1) - phi_b(i)) * exp(-tau_b(i));
  }
  
  lua_pushnumber(L, y_result);
  
  return 1;
}

/** Computes the bottom of the channel x coordinate
  
  Uses a mid-point integration rule to evaluate \f$y(\phi)\f$,
  \f[  y(\phi) = \cos(\sigma) \sum_{i=1}^{N} w_i e^{-\tau_b(\phi_{i+1/2})} \f]
  with
  
   Arguments are passed by the lua stack
 - tau_b, a table of \f$N-1\f$ values of tau at mid-points,
 - sigma, a real number \f$\sigma\f$, interior angle at base of triangle,
 - delta, a real number, the grid spacing.
 
 \param L Lua stack pointer
*/
int x_b(lua_State* L)
{
  static const int NARG = 3;
  vector tau_b = veclua_tovector(L, -NARG+0);
  vector phi_b = veclua_tovector(L, -NARG+1);
  
  real sigma = lua_tonumber(L, -NARG+2);
  
  lua_pop(L, NARG); 
  
  unsigned int n = phi_b.length();
  real x_result = 0.0;
  
  for(unsigned int i = 1; i <= (n-1); i++) {
    x_result += cos(sigma) * (phi_b(i+1) - phi_b(i)) * exp(-tau_b(i));
  }
  
  lua_pushnumber(L, x_result);
  
  return 1;
}

/** calculates x
  
  \todo needs detailed description
*/
int x_f(lua_State* L)
{
  static const int NARG = 3;
  
  vector theta = veclua_tovector(L, -NARG+0);
  vector tau   = veclua_tovector(L, -NARG+1);
  real delta   = lua_tonumber(L, -NARG+2);
  
  lua_pop(L, NARG);
     
  unsigned int n = theta.length();
  real half_length;

  real x_ip1 = 0;
  
  lua_newtable(L);
  
  lua_pushnumber(L, n);
  lua_pushnumber(L, x_ip1);
  lua_rawset(L, -3);
  
  for(unsigned int i = n-1; i >= 1; i--) {
    real x_i = x_ip1 - delta * (exp(-tau(i)) *
                                cos((theta(i) + theta(i+1)) / 2.0));
    lua_pushnumber(L, i);
    lua_pushnumber(L, x_i);
    lua_rawset(L, -3);
    
    x_ip1 = x_i;
  }

  half_length = -x_ip1 / 2.0;

  
  for(unsigned int i = 1; i <= n; i++) {
    lua_pushnumber(L, i);
    lua_rawget(L, -2);
    real x_i = lua_tonumber(L, -1);
    lua_pop(L, 1);
    x_i += half_length;
    lua_pushnumber(L, i);
    lua_pushnumber(L, x_i);
    lua_rawset(L, -3);
  }

  return 1;
}


