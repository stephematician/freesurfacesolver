#include "integrals_sbox_formulation.h"

typedef real (*_theta_func)(const real, const real_vector&);

real gaussian_theta(const real phi, const real_vector& topography) {

    const real A = topography(1);
    const real B = topography(2);
    const real L = topography(3);

    const real dy_b = -A * L * (2.0*B*B*phi) * exp(-(B*phi)*(B*phi)) / sqrt(pi);

    return atan(dy_b);

}


real smoothbox_theta(const real phi, const real_vector& topography) {

    const real A     = topography(1);
    const real B     = topography(2);
    const real L     = topography(3);
    const real phi_c = topography(4);

    const real dy_b = (A*B/(2.0*L)) *
                          (1.0 / (cosh(B*(phi+phi_c))*cosh(B*(phi+phi_c))) - 
                               1.0 / (cosh(B*(phi-phi_c))*cosh(B*(phi-phi_c))));

    return atan(dy_b);

}


void compute_convolution_s(const int                i,
                           const real_vector& theta_s,
                           const real_vector&     phi,
                           const real_vector&    dphi,
                           const real_vector& phi_mid,
                           real_vector&             o) {

    for(int j = 1; j <= theta_s.length(); j++) {
        o(j)  = (theta_s(j) - ((theta_s(i) + theta_s(i+1))/2.0)) * 
                    exp(pi*phi(j));
        o(j) /= exp(pi*phi(j)) - exp(pi*phi_mid(i));
        o(j) *= dphi(j);
    }
}


void compute_convolution_symmetric_s(const int                i,
                                     const real_vector& theta_s,
                                     const real_vector&     phi,
                                     const real_vector&    dphi,
                                     const real_vector& phi_mid,
                                     real_vector&             o) {

    for(int j = 1; j <= theta_s.length(); j++) {
        o(j)  = 0.5 * (theta_s(j) - ((theta_s(i) + theta_s(i+1))/2.0)) / 
                    tanh(pi*(phi(j) - phi_mid(i)) / 2.0);
        o(j) += 0.5 * (theta_s(j) + ((theta_s(i) + theta_s(i+1))/2.0)) / 
                    tanh(pi*(phi(j) + phi_mid(i)) / 2.0);
        o(j) *= dphi(j);
    }

}


real far_downstream_integrand_s(real x, void *_params) {

    const real * const params = (real *)_params;

    const real D         = params[0];
    const real lambda    = params[1];
    const real phi       = params[2];
    const real phi_match = params[3];
    const real gamma     = params[4];

    real f = 0.0;

    if(x != 0.0) {
        const real alpha_hat  = exp(pi*phi_match) +
                                    (1 - pow(x, gamma)) / (pow(x, gamma));
        // taken -ve to change direction of integral
        const real dalpha_hat = gamma / (pow(x, gamma+1.0)); 
        const real alpha      = -exp(pi*phi);
    
        f = atan((D*sin(lambda)*pow(alpha_hat,-lambda/pi)) /
                     (D*cos(lambda)*pow(alpha_hat,-lambda/pi)+1.0));
        f *= dalpha_hat;
        f /= (alpha_hat + alpha) * pi;
    }

    return f;

}


real far_upstream_integrand_s(real x, void *_params) {

    const real * const params = (real *)_params;

    const real D      = params[0];
    const real lambda = params[1];
    const real phi    = params[2];

    real f = 0.0;

    if(x != 0.0) {
        const real alpha_hat  = x;
        const real dalpha_hat = 1.0;
        const real alpha      = -exp(pi*phi);

        f = atan(-(D*sin(lambda)*pow(alpha_hat,lambda/pi)) /
                    (D*cos(lambda)*pow(alpha_hat,lambda/pi)+1.0));
        f *= dalpha_hat;
        f /= (alpha_hat + alpha) * pi;
    }

    return f;

}


void removed_singularity_general_s(const real_vector& theta_s,
                                   const real_vector&     phi,
                                   const real_vector& phi_mid,
                                   real_vector&             o) {

    const unsigned int n = theta_s.length();

    for(unsigned int i = 1; i <= (n-1); i++) {
        o(i) = (theta_s(i) + theta_s(i+1)) *
                   log(fabs((exp(pi*phi(n)) - exp(pi*phi_mid(i))) /
	                        (exp(pi*phi(1)) - exp(pi*phi_mid(i))))) /
                   (2.0 * pi);
  }

}


void removed_singularity_symmetric_s(const real_vector& theta_s,
                                     const real_vector&     phi,
                                     const real_vector& phi_mid,
                                     real_vector&             o) {

    const unsigned int n = theta_s.length();

    for(unsigned int i = 1; i <= (n-1); i++) {
        o(i) = (theta_s(i) + theta_s(i+1)) *
                   log(fabs(sinh(pi*(phi(n) - phi_mid(i))/2.0) / 
                               sinh(pi*(phi(n) + phi_mid(i))/2.0))) /
                   (2.0 * pi);
    }

}


real gaussian_b(real x, void *_params) {

    const real_vector topography(3, (real*)_params);

    const real B = ((real *)_params)[1];
    const real phi = ((real *)_params)[3];

    real f = 0.0;

    if(x != 0.0) {

        const real phi_hat  = atanh(x) / (B * sqrt(2) * atanh(0.5));
        const real dphi_hat = 1.0 / (B * sqrt(2) * (1.0 - x*x) * atanh(0.5));

        const real theta_b = gaussian_theta(phi_hat, topography);

        f = -theta_b * exp(pi*phi_hat);
        f *= dphi_hat;
        if (!isinf(exp(pi*phi_hat))) {
            f /= exp(pi*phi_hat) + exp(pi*phi);
        } else {
            f = 0.0;
        }
    }

    return f;

}


real smoothbox_b(real x, void *_params) {

    const real_vector topography = real_vector(4, (real*)_params);

    const real phi_c =  ((real *)_params)[3];
    const real phi   = ((real *)_params)[4];

    real f = 0.0;

    if(x != 0.0) {
        const real phi_hat  = phi_c * atanh(x) / atanh(0.5);
        const real dphi_hat = phi_c / ((1.0 - x*x) * atanh(0.5));

        const real theta_b = smoothbox_theta(phi_hat, topography);

        f = -theta_b * exp(pi*phi_hat);
        f *= dphi_hat;
        f /= exp(pi*phi_hat) + exp(pi*phi);
    }

    return !isnan(f) ? f : 0.0;

}


void topography_compute_convolution_s(const real           phi_b,
                                      const real_vector& theta_s,
                                      const real_vector&     phi,
                                      const real_vector&    dphi,
                                      real_vector&             o) {

    for(int j = 1; j <= theta_s.length(); j++) {
        o(j)  = theta_s(j) * exp(pi*phi(j));
        o(j) /= exp(pi*phi(j)) + exp(pi*phi_b);
        o(j) *= dphi(j);
    }

}


void topography_compute_symmetric_convolution_s(const real           phi_b,
                                                const real_vector& theta_s,
                                                const real_vector&     phi,
                                                const real_vector&    dphi,
                                                real_vector&             o) {

    for(int j = 1; j <= theta_s.length(); j++) {
        o(j)  = theta_s(j) * 
                    (exp(pi*phi(j))  / (exp(pi*phi(j)) + exp(pi*phi_b)) - \
                         exp(-pi*phi(j)) / (exp(-pi*phi(j)) + exp(pi*phi_b)));
        o(j) *= dphi(j);
    }

}


/* TODO: should document properly - see Chapter 4 of thesis. */
real topography_far_upstream_s(real x, void *_params) {

    const real * const params = (real *)_params;
    const real D      = params[0];
    const real lambda  = params[1];
    const real phi     = params[2];

    real f = 0.0;
    if(x != 0.0) {
        const real alpha_hat  = x;
        const real dalpha_hat = 1.0;
        const real alpha      = exp(pi*phi);
    
        f = atan(-(D*sin(lambda)*pow(alpha_hat,lambda/pi)) /
                (D*cos(lambda)*pow(alpha_hat,lambda/pi)+1.0));
        f *= dalpha_hat;
        f /= (alpha_hat + alpha) * pi;
    }

    return f;

}


/* TODO: should document properly, see Chapter 4 of thesis  */
real topography_far_downstream_s(real x, void *_params) {

    const real * const params = (real *)_params;

    const real D         = params[0];
    const real lambda    = params[1];
    const real phi       = params[2];
    const real phi_match = params[3];
    const real gamma     = params[4];

    real f = 0.0;

    if(x != 0.0) {
         //< \todo sometimes this will lose precision?
        const real alpha_hat = exp(pi*phi_match) + 
                                   (1 - pow(x, gamma)) / (pow(x, gamma));
         // taken -ve to change direction of integral
        const real dalpha_hat = gamma / (pow(x, gamma+1.0));
        const real alpha      = exp(pi*phi);

        f = atan((D*sin(lambda)*pow(alpha_hat,-lambda/pi)) /
                 (D*cos(lambda)*pow(alpha_hat,-lambda/pi)+1.0));
        f *= dalpha_hat;
        f /= (alpha_hat + alpha) * pi;
    }

    return f;

}

inline void topography_removed_singularity_b(const real phi,
                                             const real phi_m,
                                             const real_vector& topography,
                                             _theta_func theta_f, 
                                             real& o) {

    real theta_b = theta_f(phi, topography);

    o = -theta_b * log(fabs((exp( pi*phi_m) - exp(pi*phi)) / 
                                (exp(-pi*phi_m) - exp(pi*phi)))) / pi;

}


void topography_gaussian_removed_singularity_b(const real                phi,  
                                               const real              phi_m,
                                               const real_vector& topography,
                                               real&                       o) {

    topography_removed_singularity_b(       phi,          phi_m,
                                     topography, gaussian_theta,
                                              o);

}


void topography_smoothbox_removed_singularity_b(const real                phi,  
                                                const real              phi_m,
                                                const real_vector& topography,
                                                real&                       o) {

    topography_removed_singularity_b(       phi,           phi_m,
                                     topography, smoothbox_theta,
                                              o);

}


/* TODO: should document this properly - from Chapter 4 of Thesis */
inline real topography_integrand_b(const real            phi_hat,
                                   const real                phi,
                                   const real_vector& topography,
                                   _theta_func           theta_f) {

    real f = 0.0;

    if(fabs(phi_hat - phi) > 1e-9) {
        const real theta_b_hat = theta_f(phi_hat, topography);
        const real theta_b     = theta_f(    phi, topography);

        f = -(theta_b_hat - theta_b) * exp(pi*phi_hat);
        f /= exp(pi*phi_hat) - exp(pi*phi);
    } else {
        static const real mini_eps = 1e-5;
        real f1, f2;

        const real theta_b_hat1 = theta_f(phi_hat+mini_eps,topography);
        const real theta_b_hat2 = theta_f(phi_hat-mini_eps,topography);
        const real theta_b      = theta_f(             phi, topography);

        f1 = -(theta_b_hat1 - theta_b) * exp(pi*(phi_hat+mini_eps));
        f1 /= exp(pi*(phi_hat+mini_eps)) - exp(pi*phi);

        f2 = -(theta_b_hat2 - theta_b) * exp(pi*(phi_hat-mini_eps));
        f2 /= exp(pi*(phi_hat-mini_eps)) - exp(pi*phi);
        f = 0.5 * (f1 + f2);
    }

    return f;

}


real topography_gaussian_integrand_b(real x, void *_params) {

    const real_vector topography(3, (real*)_params);
    const real phi   = ((real *)_params)[3];

    const real phi_hat  = x;

    return topography_integrand_b(phi_hat, phi, topography, gaussian_theta);

}


real topography_smoothbox_integrand_b(real x, void *_params) {

    const real_vector topography(4, (real*)_params);
    const real phi   = ((real *)_params)[4];

    const real phi_hat  = x;

    return topography_integrand_b(phi_hat, phi, topography, smoothbox_theta);

}

