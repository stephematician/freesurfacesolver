/** \file integrals_sbox.cpp

  C implementation of integral formulas (e.g. trapezoidal evalution of Cauchy
  P.V. integrals). Designed for smoothed topographical disturbances.

  \todo Actually document integrals_sbox.cpp
*/

#include "integrals_sbox.h"

#include <memory>

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

    /** \todo Fix which name to use for the gsl_integration_workspace registry
              'key' */
    lua_pushstring(L, integrals_lualib_regkey);
    /* Reference to integration workspaces for this module are kept as 
     * Lua userdata */ 
    gsl_integration_workspace **ws =
        (gsl_integration_workspace**)lua_newuserdata(
                                         L,
                                         2 * sizeof(gsl_integration_workspace **)
                                     );
    ws[0] = gsl_integration_workspace_alloc(1000);
    ws[1] = gsl_integration_workspace_alloc(1000);

    /** \todo Make a new error handler. */
    gsl_set_error_handler(&integrals_sbox_handler); 

    /* Metatable for integral module - ensures clean up of integration 
     * workspaces */
    lua_newtable(L);
    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, integrals_sbox_clean_up);
    lua_settable(L, -3);
    lua_setmetatable(L, -2);
  
    lua_settable(L, LUA_REGISTRYINDEX);
  
    return 1;

}

/* Clean up function for the gsl_integration_workspace */
int integrals_sbox_clean_up(lua_State* L) {

    lua_pushstring(L, integrals_lualib_regkey);
    lua_gettable(L, LUA_REGISTRYINDEX);
    gsl_integration_workspace **ws = 
        (gsl_integration_workspace**)lua_touserdata(L, -1);
  
    gsl_integration_workspace_free(ws[0]);
    gsl_integration_workspace_free(ws[1]);

    return 0;

}

/* Two types of error code are reported to standard output then ignored 
 * all other errors result in an exit */
void integrals_sbox_handler(const char *reason, 
                            const char *file, 
                            int line, 
                            int gsl_errno) {

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

inline real average_mid_trapezoidal(const real_vector& beta,
                                    const real_vector& integrand) {

    const real n = beta.length();
    real I = 0.0f;

    for(unsigned int j = 1; j <= n; j++) {
        I += (j == 1 ?
                  beta(2) - beta(1) :
                  (j == n ?
                       beta(n) - beta(n-1) :
                       beta(j+1) - beta(j-1)
                  )
             ) * integrand(j) / 2.0;
    }

    return I;

}

void tau_s_worker(const real_vector&        theta_s,
                  const real_vector&        phi_sub,
                  const real_vector&       dphi_sub,
                  const real_vector&           beta,
                  const real_vector&       beta_sub,
                  const real_vector&     topography,
                  const ff_params&      ffus_params,
                  const ff_params&      ffds_params,
                  const real                      s,
                  gsl_integration_workspace**    ws,
                  const method_s_funcs&      method,
                  real_vector&                tau_s) {

    const unsigned int n = beta.length();

    /* Heap for vectors needed
     * name      size
     * phi        (n)
     * dphi       (n)
     * integrand  (n)
     * beta_mid (n-1)
     * phi_mid  (n-1)
     * removed  (n-1) */
    const std::unique_ptr<real[]> heap_v(new real[(3*n)+(4*(n-1))]);

    real_vector phi(n, &(heap_v[0]));
    real_vector dphi(n, &(heap_v[n]));
    real_vector integrand(n, &(heap_v[2*n]));

    real_vector beta_mid(n-1, &(heap_v[3*n]));
    real_vector phi_mid(n-1, &(heap_v[(3*n)+(n-1)]));
    real_vector removed(n-1, &(heap_v[(3*n)+(2*(n-1))]));
  
    for(unsigned int i = 1; i <= n-1; i++)
        beta_mid(i) = (beta(i) + beta(i+1)) / 2.0;

    method.transform_phi( beta, beta_sub, phi_sub, dphi_sub,  phi, s);
    method.transform_dphi(beta, beta_sub, phi_sub, dphi_sub, dphi, s);

    method.transform_phi(beta_mid, beta_sub, phi_sub, dphi_sub, phi_mid, s);

    // Set up GSL integration
 
    real downstream_p [] = {ffds_params.D,
                            ffds_params.lambda,
                            0,
                            ffds_params.phi_match,
                            ffds_params.gamma};
    real upstream_p [] = {ffus_params.D,
                          ffus_params.lambda,
                          0};

    // convert topography params into an array
    const std::unique_ptr<real[]> topography_p(new real[topography.length()+1]);
    for(int k = 1; k <= topography.length(); k++)
        topography_p[k-1] = topography(k);

    const real us_lim = exp(pi*ffus_params.phi_match);

    gsl_function downstream_F;
    gsl_function upstream_F;
    gsl_function topography_F;

    downstream_F.function = method.linear_downstream_int; 
    downstream_F.params   = downstream_p;

    upstream_F.function = method.linear_upstream_int;
    upstream_F.params   = upstream_p;  

    topography_F.function = method.topography_integrand;
    topography_F.params   = topography_p.get();

    method.removed_singularity(theta_s, phi, phi_mid, removed);

    for(unsigned int i = 1; i <= n-1; i++) {
        tau_s(i) = 0.0;

        method.convolution_terms(i, theta_s, phi, dphi, phi_mid, integrand);
 
        // Implementation of a 'trapezoidal' type integral
        tau_s(i) += average_mid_trapezoidal(beta, integrand);

        // Removal of singularity
        tau_s(i) += removed(i);

        real abserr;
        // Integration over the topography of the channel 
        real topography_term;

        topography_p[topography.length()] = phi_mid(i);
        gsl_integration_qag(&topography_F, -1.0, 1.0,      // function and limits
                            0.0, 1e-7,                     // use relative err < 1e-7
                            1000,                          // workspace size
                            GSL_INTEG_GAUSS21,             // 21 point Gauss rule
                            ws[0],                         // workspace  address
                            &topography_term, &abserr);     // return values
        tau_s(i) += topography_term;

        // Integration in the far-downstream domain 
        real ds_term;
        downstream_p[2] = phi_mid(i);
        gsl_integration_qag(&downstream_F, 0.0, 1.0,  // function and limits
                            0.0, 1e-7,                // use relative err < 1e-7
                            1000,                     // workspace size
                            GSL_INTEG_GAUSS21,        // 21 point Gauss rule
                            ws[0],                    // workspace  address
                            &ds_term, &abserr);       // return values
        tau_s(i) += ds_term;

        // Integration in the far-upstream domain 
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

}

void z_s_worker(const real_vector&  theta_mid,
                const real_vector&    tau_mid,
                const real_vector&       beta,
                const real_vector&   beta_sub,
                const real_vector&    phi_sub,
                const real_vector&   dphi_sub,
                const real                  s,
                transform_func transform_dphi,
                const cosine_or_sine  cos_sin,
                real&                      z0,
                real_vector&              z_s) {

    const unsigned int n = beta.length();
   
    const std::unique_ptr<real[]> heap_v(new real[2*(n-1)]);

    real_vector beta_mid(n-1, &(heap_v[0]));
    real_vector dphi_mid(n-1, &(heap_v[n-1]));

    for(unsigned int i = 1; i <= n-1; i++)
        beta_mid(i) = (beta(i) + beta(i+1)) / 2.0;
  
    transform_dphi(beta_mid, beta_sub, phi_sub, dphi_sub, dphi_mid, s);
  
    real M_j = 0.0;

    for(unsigned int j = n-1; j >= 1; j--) {
        const real f_j = -dphi_mid(j) *
                              exp(-tau_mid(j)) * 
                              (cos_sin == USE_COSINE ? 
                                   cos(theta_mid(j)) :
                                   sin(theta_mid(j))
                              );

        M_j    = M_j + f_j * (beta(j+1) - beta(j));
        z_s(j) = M_j - f_j * (beta(j+1) - beta(j)) / 2.0;
    }
  
    z0 = M_j;
  
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


//  for(unsigned int i = 1; i <= n; i++) {
//    real abserr;
//    
//    int int_x_b_status = gsl_integration_qag(&x_integrand_F,    // function
//                                             0, phi(i),         //  limits
//                                             0.0, 1e-7,         // relerr<1e-7
//                                             1000,              // ws size
//                                             GSL_INTEG_GAUSS21, // 21 pt Gauss
//                                             ws[0],             // workspace
//                                             &(x_b(i)),         // integral val
//                                             &abserr);          // error
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
//  }
//
//  return x_b;
//
//}

/** GSL integration to evalaute \f$y_b(\phi)\f$ at specified values of \f$\phi\f$.

  \param phi Vector of required geometry data.
  \param topography The functions used to compute the topography.
  \param theta_mid Value of \f$\theta_S\f$ at mid-points of polynomial grid.
  \param beta_sub Values of \f$\beta^{[i]}\f$ to specify polynomial grid.
  \param phi_sub Values of \f$\phi^{[i]}\f$ to specify polynomial grid.
  \param dphi_sub Values of \f$\mathrm{d}\phi^{[i]}\f$ to specify polynomial grid.
  \param s Homotopy parameter for polynomial/linear grid.

*/


void z_b_worker(const real_vector&            phi,
                const real_vector&        theta_s,
                const real_vector&         beta_s,
                const real_vector&       beta_sub,
                const real_vector&        phi_sub,
                const real_vector&       dphi_sub,
                const real                      s,
                gsl_integration_workspace**    ws,
                method_b_funcs             method,
                const real_vector&     topography,
                const ff_params&      ffus_params,
                const ff_params&      ffds_params,
                const cosine_or_sine      cos_sin,
                real_vector&                  z_b) {

    const unsigned int n = phi.length();

    gsl_function integrand_F;

    /* Should be safe with either shallow/deep copy constructor for vectors */
    z_b_integrand_params params = {/* method_b */         method,
                                   /* theta_s */         theta_s,
                                   /* beta_s */           beta_s,
                                   /* beta_sub */       beta_sub,
                                   /* phi_sub */         phi_sub,
                                   /* dphi_sub */       dphi_sub,
                                   /* ffus_params */ ffus_params,
                                   /* ffds_params */ ffds_params,
                                   /* homotopy_s */            s,
                                   /* topography */   topography,
                                   /* cos_sin */         cos_sin,
                                   /* ws */                   ws};

    integrand_F.function = &z_b_integrand;
    integrand_F.params = &params;

    for(unsigned int i = 1; i <= n; i++) {
        real abserr;
        gsl_integration_qag(&integrand_F, 0, phi(i), // function and limits
                            0.0, 1e-7,               // use relative err < 1e-7
                            1000,                    // workspace size
                            GSL_INTEG_GAUSS21,       // 21 point Gauss rule
                            ws[0],                   // workspace address
                            &(z_b(i)), &abserr);     // return values
    }

}

/* So this is an integrand for GSL to play with. */
real z_b_integrand(real phi_b, void* _params)
{

    z_b_integrand_params* params = (z_b_integrand_params*)_params;

    method_b_funcs& method = params->method_b;
  
    /* theta is easy to calculate */
    const real theta_b = method.theta_func(phi_b,
                                           params->topography);

    /* tau is the hard part */
    real tau_b = 0.0;

    const unsigned int n = params->beta_s.length();

    const std::unique_ptr<real[]> heap_v(new real[3*n]);

    real_vector phi_s(n, &(heap_v[0]));
    real_vector dphi_s(n, &(heap_v[n]));
    real_vector integrand_s(n, &(heap_v[2*n]));

    method.transform_phi(params->beta_s,
                         params->beta_sub,
                         params->phi_sub, 
                         params->dphi_sub, 
                         phi_s,
                         params->homotopy_s);
    method.transform_dphi(params->beta_s,
                          params->beta_sub,
                          params->phi_sub,
                          params->dphi_sub,
                          dphi_s,
                          params->homotopy_s);

    method.convolution_s(phi_b,
                         params->theta_s,
                         phi_s,
                         dphi_s,
                         integrand_s);

    tau_b += average_mid_trapezoidal(params->beta_s,
                                     integrand_s);

    gsl_function upstream_F;
    gsl_function downstream_F;

    const real us_lim = exp(pi*params->ffus_params.phi_match);

    real upstream_p[] = {params->ffus_params.D,
                         params->ffus_params.lambda,
                         phi_b};
    real downstream_p[] = {params->ffds_params.D,
                           params->ffds_params.lambda,
                           phi_b,
                           params->ffds_params.phi_match,
                           params->ffds_params.gamma};
  
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
                        params->ws[1],               // workspace size address
                        &us_term, &abserr);          // return values
    tau_b += us_term;
  
    real ds_term;
    gsl_integration_qag(&downstream_F, 0.0, 1.0,     // function and limits
                        0.0, 1e-7,                   // use relative err < 1e-7
                        1000,                        // workspace size
                        GSL_INTEG_GAUSS21,           // 21 point Gauss rule
                        params->ws[1],               // workspace size address
                        &ds_term, &abserr);          // return values
    tau_b += ds_term;
  
    /* Convert topography params into an array */
    const std::unique_ptr<real[]> topography_p(
                                      new real[params->topography.length()+1]
                                  );
    for(int k = 1; k <= params->topography.length(); k++)
        topography_p[k-1] = params->topography(k);
    topography_p[params->topography.length()] = phi_b;

    gsl_function integrand_b_F;

    integrand_b_F.function = method.topography_integrand;
    integrand_b_F.params   = topography_p.get();

    real phi_m = fmax(fabs(phi_s(1)), fabs(phi_s(n)));
    real phi_M = fmax(phi_m, fabs(phi_b*2.0)); // to avoid issues with singularity

    abserr = 0;
    real integral_b = 0.0;
    int int_b_status = gsl_integration_qag(&integrand_b_F,    // function
                                           -phi_M, phi_M,     // limits
                                           0.0, 1e-8,         // rel err < 1e-8
                                           1000,              // workspace size
                                           GSL_INTEG_GAUSS21, // 21 pt Gauss rule
                                           params->ws[1],     // workspace
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

    real removed_singularity_b = 0.0;
    method.removed_singularity_b(phi_b,
                                 phi_M,
                                 params->topography,
                                 removed_singularity_b);
    tau_b += removed_singularity_b;

    return exp(-tau_b) * (params->cos_sin == USE_COSINE ?
                              cos(theta_b) :
                              sin(theta_b));

}

