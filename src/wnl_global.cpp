// #define LUA_LIB
extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "vector.h"
#include "gsl/gsl_odeiv.h"

#define pi M_PI

extern "C" int luaopen_wnl_global(lua_State *);
static int wnl_compute(lua_State*);
int wnl_dx(double, const double *, double *, void *); 
int wnl_J(double, const double *, double *, double *, void *);

static const luaL_reg wnlR[] = {
{"compute", wnl_compute},
{NULL, NULL}
};

static int wnl_compute(lua_State* L)
{

  typedef struct {
    double A;
    double y_inf;
    double FROUDE;
  } params_struct;
  
  union {
    double raw[3];
    params_struct s;
  } params;
  
  double y[2], yerr[2], dydx_in[2], dydx_out[2];
  
  params.s.A = lua_tonumber(L, -5);
  params.s.y_inf = lua_tonumber(L, -4);
  params.s.FROUDE = lua_tonumber(L, -3);
  
  double N = lua_tonumber(L, -2);
  double delta = lua_tonumber(L, -1);
  
  lua_pop(L, 5);

  vector theta((int)N); // = vector((int)N);
   
  gsl_odeiv_system wnl_sys;
  
  wnl_sys.function = &wnl_dx;
  wnl_sys.jacobian = &wnl_J;
  wnl_sys.dimension = 2;
  wnl_sys.params = params.raw;
  
  gsl_odeiv_step *ode45 = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, 2);
  
  double x = -delta * ((N - 1.0) / 2.0);
  y[0] = params.s.y_inf - 1;
  y[1] = 0.0;
  dydx_in[0] = 0.0;
  dydx_in[1] = 0.0;
  
  for(unsigned int i = 2; i <= N; i++) {
    theta(i-1) = atan(y[1]);
    int rv = gsl_odeiv_step_apply(ode45, x, delta, y, yerr, dydx_in, dydx_out, &wnl_sys);
    dydx_in[0] = dydx_out[0];
    dydx_in[1] = dydx_out[1];
    x += delta;
  }
  theta(N) = atan(y[1]);

  veclua_pushtable(L, theta);

  gsl_odeiv_step_free(ode45);
  
  return 1;
}

extern "C" int luaopen_wnl_global(lua_State *L) {
  luaL_register(L, "wnl", wnlR);
  return 1;
}

double delta_approx(double x, double B)
{
  return B * exp(-B * B * x * x) / sqrt(pi);
}

int wnl_dx(double x, const double *y, double *dydx, void *params) 
{
  double eta = y[0];
  double zeta = y[1];
  double A = ((double*)params)[0];
  double F = ((double*)params)[2];
  
  dydx[0] = zeta;
  dydx[1] = eta * (6.0 * (F - 1.0) + (-9.0 * eta / 2.0)) + 
            (-3.0 * A * ((x < 0.0) ? 1.0 : 0.0));
  
  return GSL_SUCCESS; 
}

int wnl_J(double x, const double *y, double *dfdy, double *dfdx, void *params) 
{
  double eta = y[0];
  double zeta = y[1];
  double A = ((double*)params)[0];
  double yinf = ((double*)params)[1];
  double F = ((double*)params)[2];
  
  dfdx[0] = 0.0;
  dfdx[1] = -(3.0 * A * delta_approx(x, 2.0));
  
  dfdy[0] = 0.0;
  dfdy[1] = 1.0;
  dfdy[2] = 6.0 * (F - 1.0) - (9.0 * eta);
  dfdy[3] = 0.0;
  
  return GSL_SUCCESS; 
}

