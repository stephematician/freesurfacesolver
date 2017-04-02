/** \file integrals_sbox.cpp

  C implementation of integral formulas (e.g. trapezoidal evalution of Cauchy
  P.V. integrals). Designed for smoothed topographical disutrbances.

  \todo Actually document integrals_sbox.cpp
*/

#include "integrals_sbox.h"

static int integrals_sbox_clean_up(lua_State *);

extern char integrals_lualib_name[];
extern char integrals_lualib_regkey[];

extern luaL_Reg integralsR[];

int luaopen_integrals_sbox(lua_State *L) {
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
  gsl_integration_workspace **ws = (gsl_integration_workspace**)lua_newuserdata(L, 2*sizeof(gsl_integration_workspace **));
  ws[0] = gsl_integration_workspace_alloc(1000);
  ws[1] = gsl_integration_workspace_alloc(1000);

  /** \todo Make a new error handler. */
  gsl_set_error_handler(&integrals_sbox_handler); 

  // Make a new metatable
  lua_newtable(L);
  // Garbage collection method
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, integrals_sbox_clean_up);
  lua_settable(L, -3);
  // Set metatable
  lua_setmetatable(L, -2);
  
  // Should now have stack (T->B) -> userdata, "integrals_sbox", luaL register table
  lua_settable(L, LUA_REGISTRYINDEX);
  
 // Should now have stack (T->B) -> luaL register table  
  return 1;
}

/* Clean up function for the gsl_integration_workspace */
int integrals_sbox_clean_up(lua_State* L)
{
  lua_pushstring(L, integrals_lualib_regkey);
  lua_gettable(L, LUA_REGISTRYINDEX);
  gsl_integration_workspace **ws = (gsl_integration_workspace**)lua_touserdata(L, -1);
  
  gsl_integration_workspace_free(ws[0]);
  gsl_integration_workspace_free(ws[1]);

  return 0;
}

void integrals_sbox_handler(const char *reason, 
                            const char *file, 
                            int line, 
                            int gsl_errno)
{
  switch(gsl_errno) {
    case 0 :
    // This shouldn't happen, right?
    break;
    case GSL_EROUND :
      printf("Warning: GSL was not able to integrate to desired accuracy.\n");
      printf("%s\n", reason);
      printf("file %s\n", file);
      printf("line %u\n", line);
    break;
    case GSL_EMAXITER :
      printf("Warning: GSL was not able to integrate to desired accuracy.\n");
      printf("%s\n", reason);
      printf("file %s\n", file);
      printf("line %u\n", line);
    break;
    default:
      exit(EXIT_FAILURE);
  }
}

vector tau_s_worker(const vector&             theta_s,
                    const vector&             phi_sub,
                    const vector&            dphi_sub,
                    const vector&                beta,
                    const vector&            beta_sub,
                    const vector&          topography,
                    const linear_params&  lin_ustream,
                    const linear_params&  lin_dstream,
                    const real                      s,
                    gsl_integration_workspace**    ws,
                    method_s_funcs&            method)
{

  unsigned int n = beta.length();

  vector phi  = vector(n);
  vector dphi = vector(n);

  vector integrand = vector(n);

  vector beta_mid = vector(n-1);
  vector phi_mid  = vector(n-1);

  vector tau_s = vector(n-1);
  
  for(unsigned int i = 1; i <= n-1; i++) {
    beta_mid(i) = (beta(i) + beta(i+1)) / 2.0;
  }

  (*(method.transform_phi))(beta, beta_sub, phi_sub, dphi_sub, phi, s);
  (*(method.transform_dphi))(beta, beta_sub, phi_sub, dphi_sub, dphi, s);
  
  (*(method.transform_phi))(beta_mid, beta_sub, phi_sub, dphi_sub, phi_mid, s);
  
  // Set up GSL integration
  
  real* downstream_p = new real[5];
  real* upstream_p = new real[3];
  real* topography_p = new real[topography.length()+1];
  
  upstream_p[0] = lin_ustream.D;      downstream_p[0] = lin_dstream.D;
  upstream_p[1] = lin_ustream.lambda; downstream_p[1] = lin_dstream.lambda;
  
  downstream_p[3] = lin_dstream.phi_match;
  downstream_p[4] = lin_dstream.gamma;
  
  // convert topography params into an array
  for(int k = 1; k <= topography.length(); k++) {
    topography_p[k-1] = topography(k);
  }

  real us_lim     = exp(pi*lin_ustream.phi_match); // exp(pi*phi(1));

  gsl_function downstream_F;
  gsl_function upstream_F;
  gsl_function topography_F;

  downstream_F.function = method.linear_downstream_int; // &far_downstream_integrand_s;
  downstream_F.params = (void *)downstream_p;

  upstream_F.function = method.linear_upstream_int; // &far_upstream_integrand_s;
  upstream_F.params = (void *)upstream_p;  

  topography_F.function = method.topography_integrand;
  topography_F.params = (void *)topography_p;

  vector removed_singularity = vector(n-1);
  (*(method.removed_singularity))(theta_s, phi, phi_mid, removed_singularity);
  
  /** \todo lets get started on the topography integral code */
  
  for(unsigned int i = 1; i <= n-1; i++) {
    tau_s(i) = 0.0;

    (*(method.convolution_terms))(i, theta_s, phi, dphi, phi_mid, integrand);
  
    // This implementation of a 'trapezoidal' type integral has been tested
    for(unsigned int j = 1; j <= n; j++) {
      real delta = (j == 1 ? beta(2) - beta(1) :
                             (j == n ? beta(n) - beta(n-1) :
                                       beta(j+1) - beta(j-1)));
      tau_s(i) += delta * integrand(j) / 2.0;
    }
    // Expression for removal of singularity
    tau_s(i) += removed_singularity(i);

    real abserr;
    // Expression for integration over topography of channel
    real topography_term;

    topography_p[topography.length()] = phi_mid(i);
    gsl_integration_qag(&topography_F, -1.0, 1.0,      // function and limits
                        0.0, 1e-7,                     // use relative err < 1e-7
                        1000,                          // workspace size
                        GSL_INTEG_GAUSS21,             // 21 point Gauss rule
                        ws[0],                         // workspace  address
                        &topography_term, &abserr);     // return values
    tau_s(i) += topography_term;
  
    // Expression for integration far downstream
    real ds_term;
    downstream_p[2] = phi_mid(i);
    gsl_integration_qag(&downstream_F, 0.0, 1.0,  // function and limits
                        0.0, 1e-7,                // use relative err < 1e-7
                        1000,                     // workspace size
                        GSL_INTEG_GAUSS21,        // 21 point Gauss rule
                        ws[0],                    // workspace  address
                        &ds_term, &abserr);       // return values
    tau_s(i) += ds_term;
    
    // Expression for integration far upstream
    real us_term;
    upstream_p[2] = phi_mid(i);
    gsl_integration_qag(&upstream_F, 0, us_lim,  // function and limits
                        0.0, 1e-7,               // use relative err < 1e-7
                        1000,                    // workspace size
                        GSL_INTEG_GAUSS21,       // 21 point Gauss rule
                        ws[0],                   // workspace  address
                        &us_term, &abserr);      // return values
    tau_s(i) += us_term;
  }

  delete [] downstream_p;
  delete [] upstream_p;
  delete [] topography_p;

  return tau_s;
}

vector y_s_worker(const vector&       theta_mid,
                  const vector&         tau_mid,
                  const vector&            beta,
                  const vector&        beta_sub,
                  const vector&         phi_sub,
                  const vector&        dphi_sub,
                  const real                  s,
                  transform_func transform_dphi,
                  real&                    eta0)
{
  unsigned int n = beta.length();
    
  vector beta_mid = vector(n-1);
  for(unsigned int i = 1; i <= n-1; i++) {
    beta_mid(i) = (beta(i) + beta(i+1)) / 2.0;
  }
  
  vector dphi_mid = vector(n-1);
  (*transform_dphi)(beta_mid, beta_sub, phi_sub, dphi_sub, dphi_mid, s);
  
  vector y_s = vector(n-1);
  
  real f_j = -dphi_mid(n-1) * exp(-tau_mid(n-1)) * sin(theta_mid(n-1));
  real M_j = f_j * (beta(n) - beta(n-1));
  y_s(n-1) = M_j / 2.0;

  for(unsigned int j = n-2; j >= 1; j--) {
    f_j = -dphi_mid(j) * exp(-tau_mid(j)) * sin(theta_mid(j));
    M_j = M_j + f_j * (beta(j+1) - beta(j));
    
    y_s(j) = M_j - f_j * (beta(j+1) - beta(j)) / 2.0;
  }
  
  eta0 = M_j;
  
  return y_s;
}


vector x_s_worker(const vector&       theta_mid, 
                  const vector&         tau_mid, 
                  const vector&            beta,
                  const vector&        beta_sub,
                  const vector&         phi_sub,
                  const vector&        dphi_sub,
                  const real                  s,
                  transform_func transform_dphi,
                  real&                      x0)
{
  
  unsigned int n = beta.length();
    
  vector beta_mid = vector(n-1);
  for(unsigned int i = 1; i <= n-1; i++) {
    beta_mid(i) = (beta(i) + beta(i+1)) / 2.0;
  }
  
  vector dphi_mid = vector(n-1);
  (*transform_dphi)(beta_mid, beta_sub, phi_sub, dphi_sub, dphi_mid, s);
  
  vector x_s = vector(n-1);
  
  real f_j = -dphi_mid(n-1) * exp(-tau_mid(n-1)) * cos(theta_mid(n-1)); 
  real M_j = f_j * (beta(n) - beta(n-1));
  x_s(n-1) = M_j / 2.0;

  for(unsigned int j = n-2; j >= 1; j--) {
    f_j = -dphi_mid(j) * exp(-tau_mid(j)) * cos(theta_mid(j));
    M_j = M_j + f_j * (beta(j+1) - beta(j));
    
    x_s(j) = M_j - f_j * (beta(j+1) - beta(j)) / 2.0;
  }

  x0 = M_j;

  return x_s;
}

/** GSL integration to evalaute \f$x_b(\phi)\f$ at specified values of \f$\phi\f$.

  \param phi Vector of required geometry data.
  \param topography The functions used to compute the topography.
  \param theta_mid Value of \f$\theta_S\f$ at mid-points of polynomial grid.
  \param beta_sub Values of \f$\beta^{[i]}\f$ to specify polynomial grid.
  \param phi_sub Values of \f$\phi^{[i]}\f$ to specify polynomial grid.
  \param dphi_sub Values of \f$\mathrm{d}\phi^{[i]}\f$ to specify polynomial grid.
  \param s Homotopy parameter for polynomial/linear grid.

*/
vector x_b_worker(const vector&                 phi,
                  const vector&             theta_s,
                  const vector&              beta_s,
                  const vector&            beta_sub,
                  const vector&             phi_sub,
                  const vector&            dphi_sub,
                  const real                      s,
                  gsl_integration_workspace**    ws,
                  method_b_funcs             method,
                  const vector&          topography,
                  const linear_params&  lin_ustream,
                  const linear_params&  lin_dstream)
{

  unsigned int n = phi.length();

  vector x_b = vector(n);
  gsl_function x_integrand_F;

  x_b_integrand_params params;

  params.method_b = method;
  params.theta_s = theta_s;
  params.beta_s = beta_s;

  params.phi_sub = phi_sub;
  params.dphi_sub = dphi_sub;
  params.beta_sub = beta_sub;

  params.ws = ws;

  params.lin_ustream = lin_ustream;
  params.lin_dstream = lin_dstream;

  params.topography = topography;
  params.homotopy_s = s;

  x_integrand_F.function = &x_b_integrand;
  x_integrand_F.params = &params;

  for(unsigned int i = 1; i <= n; i++) {
    real abserr;
    
    int int_x_b_status = gsl_integration_qag(&x_integrand_F,    // function
                                             0, phi(i),         //  limits
                                             0.0, 1e-7,         // relerr<1e-7
                                             1000,              // ws size
                                             GSL_INTEG_GAUSS21, // 21 pt Gauss
                                             ws[0],             // workspace
                                             &(x_b(i)),         // integral val
                                             &abserr);          // error
    //switch (int_x_b_status) {
    //  case 0 :
    //  break;
    //  case GSL_EROUND :
    //    printf("Warning: could not integrate over channel topography to desired accuracy.\n");
    //    printf("Return value = %f\n", x_b(i));
    //    printf("Error value = %f\n", abserr);
    //  break;
    //  case GSL_EMAXITER :
    //    printf("Warning: could not integrate over channel topography to desired accuracy.\n");
    //    printf("Return value = %f\n", x_b(i));
    //    printf("Error value = %f\n", abserr);
    //  break;
    //  default:
    //    return NAN;
    //}
  }

  return x_b;

}

/* So this is an integrand for GSL to play with. Hopefully x_b_worker is all good. */
real x_b_integrand(real phi_b, void* params)
{
  method_b_funcs method = (*((x_b_integrand_params*)params)).method_b; // ???
  
  vector theta_s = (*((x_b_integrand_params*)params)).theta_s;  // ???
  vector beta_s = (*((x_b_integrand_params*)params)).beta_s;  // ???

  vector beta_sub = (*((x_b_integrand_params*)params)).beta_sub;
  vector phi_sub = (*((x_b_integrand_params*)params)).phi_sub;
  vector dphi_sub = (*((x_b_integrand_params*)params)).dphi_sub;

  real s = (*((x_b_integrand_params*)params)).homotopy_s;

  vector topography = (*((x_b_integrand_params*)params)).topography; //  ???

  gsl_integration_workspace** ws = (*((x_b_integrand_params*)params)).ws; // ???

  /* theta is easy to calculate */
  real theta_b = (*(method.theta_func))(phi_b, topography);

  /* tau is the hard part */
  real tau_b = 0.0;

  /** \todo Check contribution due to free-surface on finite interval \f$\phi^{[0]} \to \phi^{[M]}\f$. */
  unsigned int n = beta_s.length();

  vector phi_s  = vector(n);
  vector dphi_s = vector(n);

  vector integrand_s = vector(n);

  /*vector tau_s = vector(n-1);*/

  /*for(unsigned int i = 1; i <= n-1; i++) {
    beta_mid(i) = (beta(i) + beta(i+1)) / 2.0;
    }
  */

  (*(method.transform_phi))(beta_s, beta_sub, phi_sub, dphi_sub, phi_s, s);
  (*(method.transform_dphi))(beta_s, beta_sub, phi_sub, dphi_sub, dphi_s, s);

  /* I need to check if this is ok? i may want to change the first term? */
  (*(method.convolution_s))(phi_b, theta_s, phi_s, dphi_s, integrand_s);

  /* This implementation of a 'trapezoidal' type integral has been tested */
  for(unsigned int j = 1; j <= n; j++) {
    real delta = (j == 1 ? beta_s(2) - beta_s(1) :
                           (j == n ? beta_s(n) - beta_s(n-1) :
                                     beta_s(j+1) - beta_s(j-1)));
    tau_b += delta * integrand_s(j) / 2.0;
  }

  /** \todo Check contribution due to free-surface due to linearised downstream approximation. */
  /** \todo Check contribution due to free-surface due to linearised upstream approximation. */
  
  gsl_function upstream_F;
  gsl_function downstream_F;

  real* upstream_p = new real[3];
  real* downstream_p = new real[5];

  linear_params lin_ustream = (*((x_b_integrand_params*)params)).lin_ustream;
  linear_params lin_dstream = (*((x_b_integrand_params*)params)).lin_dstream;

  real us_lim = exp(pi*lin_ustream.phi_match);

  upstream_p[0] = lin_ustream.D;      downstream_p[0] = lin_dstream.D;
  upstream_p[1] = lin_ustream.lambda; downstream_p[1] = lin_dstream.lambda;
  upstream_p[2] = phi_b;              downstream_p[2] = phi_b;

  downstream_p[3] = lin_dstream.phi_match;
  downstream_p[4] = lin_dstream.gamma;

  upstream_F.function = method.linear_upstream_s;
  upstream_F.params = upstream_p;

  downstream_F.function = method.linear_downstream_s;
  downstream_F.params = downstream_p;

  real abserr = 0.0;
  real us_term;
  gsl_integration_qag(&upstream_F, 0.0, us_lim,    // function and limits
                      0.0, 1e-7,                   // use relative err < 1e-7
                      1000,                        // workspace size
                      GSL_INTEG_GAUSS21,           // 21 point Gauss rule
                      ws[1],                       // workspace size address
                      &us_term, &abserr);          // return values
  tau_b += us_term;
  
  real ds_term;
  gsl_integration_qag(&downstream_F, 0.0, 1.0,     // function and limits
                      0.0, 1e-7,                   // use relative err < 1e-7
                      1000,                        // workspace size
                      GSL_INTEG_GAUSS21,           // 21 point Gauss rule
                      ws[1],                       // workspace size address
                      &ds_term, &abserr);          // return values
  //  printf("us = %.5e, ds = %.5e\n", us_term, ds_term);
  tau_b += ds_term;
  
  /** \todo Check this main contribution due to convolution on bottom */
  real* topography_p = new real[topography.length()+1];

  /* Convert topography params into an array */
  for(int k = 1; k <= topography.length(); k++) {
    topography_p[k-1] = topography(k);
  }
  topography_p[topography.length()] = phi_b; // I hope.

  gsl_function integrand_b_F;

  integrand_b_F.function = method.topography_integrand;
  integrand_b_F.params   = topography_p;

  real phi_m = fmax(fabs(phi_s(1)), fabs(phi_s(n)));
  real phi_M = fmax(phi_m, fabs(phi_b*2.0)); // to avoid issues with singularity

  real integral_b = 0.0;

  int int_b_status = gsl_integration_qag(&integrand_b_F,    // function
                                        -phi_M, phi_M,      // limits
                                         0.0, 1e-8,         // rel err < 1e-8
                                         1000,              // workspace size
                                         GSL_INTEG_GAUSS21, // 21 pt Gauss rule
                                         ws[1],             // workspace
                                         &integral_b,       // integral value
                                         &abserr);          // error
  switch (int_b_status) {
    case 0 :
    break;
    case GSL_EROUND :
      printf("Warning: could not integrate over channel topography to desired accuracy.\n");
      printf("Return value = %f\n", integral_b);
      printf("Error value = %f\n", abserr);
    break;
    case GSL_EMAXITER :
      printf("Warning: could not integrate over channel topography to desired accuracy.\n");
      printf("Return value = %f\n", integral_b);
      printf("Error value = %f\n", abserr);
    break;
    default:
      return NAN;
  }
  tau_b += integral_b;

  // printf("integral_b = %.5f\n", integral_b);
  /* CHECK THE SHIT OUTTA ALL OF THIS! ^^^ */

  /** \todo Check this contribution from removal of singularity. */
  real removed_singularity_b = 0.0;
  (*(method.removed_singularity_b))(phi_b,
                                    phi_M,
                                    topography,
                                    removed_singularity_b);
  tau_b += removed_singularity_b;
  // printf("singularity_b = %.5f\n", removed_singularity_b);


  real f = exp(-tau_b) * cos(theta_b);
  //  printf("theta = %.5f, tau = %.5f, integrand = %.5f\n", theta_b, tau_b, f);

  delete [] upstream_p;
  delete [] downstream_p;
  delete [] topography_p;

  return f;

}

/** GSL integration to evalaute \f$x_b(\phi)\f$ at specified values of \f$\phi\f$.

  \param phi Vector of required geometry data.
  \param topography The functions used to compute the topography.
  \param theta_mid Value of \f$\theta_S\f$ at mid-points of polynomial grid.
  \param beta_sub Values of \f$\beta^{[i]}\f$ to specify polynomial grid.
  \param phi_sub Values of \f$\phi^{[i]}\f$ to specify polynomial grid.
  \param dphi_sub Values of \f$\mathrm{d}\phi^{[i]}\f$ to specify polynomial grid.
  \param s Homotopy parameter for polynomial/linear grid.

*/
vector y_b_worker(const vector&                 phi,
                  const vector&             theta_s,
                  const vector&              beta_s,
                  const vector&            beta_sub,
                  const vector&             phi_sub,
                  const vector&            dphi_sub,
                  const real                      s,
                  gsl_integration_workspace**    ws,
                  method_b_funcs             method,
                  const vector&          topography,
                  const linear_params&  lin_ustream,
                  const linear_params&  lin_dstream)
{

  unsigned int n = phi.length();

  vector y_b = vector(n);
  gsl_function y_integrand_F;

  x_b_integrand_params params;

  params.method_b = method;
  params.theta_s = theta_s;
  params.beta_s = beta_s;

  params.phi_sub = phi_sub;
  params.dphi_sub = dphi_sub;
  params.beta_sub = beta_sub;

  params.ws = ws;

  params.lin_ustream = lin_ustream;
  params.lin_dstream = lin_dstream;

  params.topography = topography;
  params.homotopy_s = s;


  y_integrand_F.function = &y_b_integrand;
  y_integrand_F.params = &params;

  for(unsigned int i = 1; i <= n; i++) {
    real abserr;
    params.temp_phi = phi(i);
    /*if(fabs(phi(i) - 5.0) < 1e-5) printf("Test 0\n");*/
    gsl_integration_qag(&y_integrand_F, 0, phi(i), // function and limits
                        0.0, 1e-7,                 // use relative err < 1e-7
                        1000,                      // workspace size
                        GSL_INTEG_GAUSS21,         // 21 point Gauss rule
                        ws[0],                     // workspace address
                        &(y_b(i)), &abserr);       // return values
    /*if(fabs(phi(i) - 5.0) < 1e-5) printf("Test clear\n");*/
  }
  return y_b;

}

/* So this is an integrand for GSL to play with. Hopefully x_b_worker is all good. */
real y_b_integrand(real phi_b, void* params)
{
  method_b_funcs method = (*((x_b_integrand_params*)params)).method_b; // ???
  
  vector theta_s = (*((x_b_integrand_params*)params)).theta_s;  // ???
  vector beta_s = (*((x_b_integrand_params*)params)).beta_s;  // ???

  vector beta_sub = (*((x_b_integrand_params*)params)).beta_sub;
  vector phi_sub = (*((x_b_integrand_params*)params)).phi_sub;
  vector dphi_sub = (*((x_b_integrand_params*)params)).dphi_sub;

  real s = (*((x_b_integrand_params*)params)).homotopy_s;

  real temp_phi = (*((x_b_integrand_params*)params)).temp_phi;

  vector topography = (*((x_b_integrand_params*)params)).topography; //  ???

  gsl_integration_workspace** ws = (*((x_b_integrand_params*)params)).ws; // ???

  /* theta is easy to calculate */
  real theta_b = (*(method.theta_func))(phi_b, topography);

  /* tau is the hard part */
  real tau_b = 0.0;

  /** \todo Check contribution due to free-surface on finite interval \f$\phi^{[0]} \to \phi^{[M]}\f$. */
  unsigned int n = beta_s.length();

  vector phi_s  = vector(n);
  vector dphi_s = vector(n);

  vector integrand_s = vector(n);

  /*vector tau_s = vector(n-1);*/

  /*for(unsigned int i = 1; i <= n-1; i++) {
    beta_mid(i) = (beta(i) + beta(i+1)) / 2.0;
    }
  */

  (*(method.transform_phi))(beta_s, beta_sub, phi_sub, dphi_sub, phi_s, s);
  (*(method.transform_dphi))(beta_s, beta_sub, phi_sub, dphi_sub, dphi_s, s);

  /* I need to check if this is ok? i may want to change the first term? */
  (*(method.convolution_s))(phi_b, theta_s, phi_s, dphi_s, integrand_s);

  /* This implementation of a 'trapezoidal' type integral has been tested */
  for(unsigned int j = 1; j <= n; j++) {
    real delta = (j == 1 ? beta_s(2) - beta_s(1) :
                           (j == n ? beta_s(n) - beta_s(n-1) :
                                     beta_s(j+1) - beta_s(j-1)));
    tau_b += delta * integrand_s(j) / 2.0;
  }

  /** \todo Check contribution due to free-surface due to linearised downstream approximation. */
  /** \todo Check contribution due to free-surface due to linearised upstream approximation. */
  
  gsl_function upstream_F;
  gsl_function downstream_F;

  real* upstream_p = new real[3];
  real* downstream_p = new real[5];

  linear_params lin_ustream = (*((x_b_integrand_params*)params)).lin_ustream;
  linear_params lin_dstream = (*((x_b_integrand_params*)params)).lin_dstream;

  real us_lim = exp(pi*lin_ustream.phi_match);

  upstream_p[0] = lin_ustream.D;      downstream_p[0] = lin_dstream.D;
  upstream_p[1] = lin_ustream.lambda; downstream_p[1] = lin_dstream.lambda;
  upstream_p[2] = phi_b;              downstream_p[2] = phi_b;

  downstream_p[3] = lin_dstream.phi_match;
  downstream_p[4] = lin_dstream.gamma;

  upstream_F.function = method.linear_upstream_s;
  upstream_F.params = upstream_p;

  downstream_F.function = method.linear_downstream_s;
  downstream_F.params = downstream_p;

  real abserr = 0.0;
  real us_term;
  /*if(fabs(temp_phi - 5.0) < 1e-5) printf("Test 1\n");*/
  gsl_integration_qag(&upstream_F, 0.0, us_lim,    // function and limits
                      0.0, 1e-7,                   // use relative err < 1e-7
                      1000,                        // workspace size
                      GSL_INTEG_GAUSS21,           // 21 point Gauss rule
                      ws[1],                       // workspace size address
                      &us_term, &abserr);          // return values
  tau_b += us_term;
  
  real ds_term;
  /*if(fabs(temp_phi - 5.0) < 1e-5) printf("Test 2\n");*/
  gsl_integration_qag(&downstream_F, 0.0, 1.0,     // function and limits
                      0.0, 1e-7,                   // use relative err < 1e-7
                      1000,                        // workspace size
                      GSL_INTEG_GAUSS21,           // 21 point Gauss rule
                      ws[1],                       // workspace size address
                      &ds_term, &abserr);          // return values
  //  printf("us = %.5e, ds = %.5e\n", us_term, ds_term);
  tau_b += ds_term;
  
  /** \todo Check this main contribution due to convolution on bottom */
  real* topography_p = new real[topography.length()+1]; // +2 for temp_phi

  /* Convert topography params into an array */
  for(int k = 1; k <= topography.length(); k++) {
    topography_p[k-1] = topography(k);
  }
  topography_p[topography.length()] = phi_b; // I hope.
  // topography_p[topography.length()+1] = temp_phi; // TEMPORARY

  gsl_function integrand_b_F;

  integrand_b_F.function = method.topography_integrand;
  integrand_b_F.params   = topography_p;

  real phi_m = fmax(fabs(phi_s(1)), fabs(phi_s(n)));
  real phi_M = fmax(phi_m, fabs(phi_b*2.0)); // to avoid issues with singularity

  abserr = 0;
  real integral_b = 0.0;

  /*if(fabs(temp_phi - 5.0) < 1e-5) {
  printf("Test 3\n");
  printf("phi_m = %.16f\n", phi_M);
  printf("phi_M = %.16f\n", phi_M);
  }*/

  int int_b_status = gsl_integration_qag(&integrand_b_F,    // function
                                         -phi_M, phi_M,     // limits
                                         0.0, 1e-8,         // rel err < 1e-8
                                         1000,              // workspace size
                                         GSL_INTEG_GAUSS21, // 21 pt Gauss rule
                                         ws[1],             // workspace
                                         &integral_b,       // integral value
                                         &abserr);          // error
  switch (int_b_status) {
    case 0 :
    break;
    case GSL_EROUND :
      printf("Warning: could not integrate over channel topography to desired accuracy.\n");
      printf("Return value = %f\n", integral_b);
      printf("Error value = %f\n", abserr);
    break;
    case GSL_EMAXITER :
      printf("Warning: could not integrate over channel topography to desired accuracy.\n");
      printf("Return value = %f\n", integral_b);
      printf("Error value = %f\n", abserr);
    break;
    default:
      return NAN;
   }

  tau_b += integral_b;

  // printf("integral_b = %.5f\n", integral_b);
  /* CHECK THE SHIT OUTTA ALL OF THIS! ^^^ */

  /** \todo Check this contribution from removal of singularity. */
  real removed_singularity_b = 0.0;
  (*(method.removed_singularity_b))(phi_b,
                                    phi_M,
                                    topography,
                                    removed_singularity_b);
  tau_b += removed_singularity_b;
  // printf("singularity_b = %.5f\n", removed_singularity_b);

  real f = exp(-tau_b) * sin(theta_b);
  //  printf("theta = %.5f, tau = %.5f, integrand = %.5f\n", theta_b, tau_b, f);
  /*if(fabs(temp_phi - 5.0) < 1e-5) printf("Test 4 %f\n", f);
    if(fabs(temp_phi - 5.0) < 1e-5) printf("\n\n");*/

  delete []  upstream_p;
  delete [] downstream_p;
  delete [] topography_p;

  return f;

}
