#include "integrals_sbox_formulation.h"

real gaussian_theta(const real phi, const vector& topography)
{
  real A = topography(1);
  real B = topography(2);
  real L = topography(3);

  real dy_b = -A * L * (2.0*B*B*phi) * exp(-(B*phi)*(B*phi)) / sqrt(pi);

  real theta_b = atan(dy_b);

  return theta_b;
}

real smoothbox_theta(const real phi, const vector& topography)
{
  real A     = topography(1);
  real B     = topography(2);
  real L     = topography(3);
  real phi_c = topography(4);

  real dy_b = (A*B/(2.0*L)) * (1.0 / \
                               (cosh(B*(phi+phi_c))*cosh(B*(phi+phi_c))) - \
                               1.0 / \
                               (cosh(B*(phi-phi_c))*cosh(B*(phi-phi_c))));
  real theta_b = atan(dy_b);

  return theta_b;
}

void compute_convolution_s(const int           i,
                           const vector& theta_s,
                           const vector&     phi,
                           const vector&    dphi,
                           const vector& phi_mid,
                           vector&             o)
{
  for(int j = 1; j <= theta_s.length(); j++) {
    o(j)  = (theta_s(j) - ((theta_s(i) + theta_s(i+1))/2.0)) * exp(pi*phi(j));
    o(j) /= exp(pi*phi(j)) - exp(pi*phi_mid(i));
    o(j) *= dphi(j);
  }
}

real far_downstream_integrand_s(real x, void *params)
{
  real D         = ((real *)params)[0];
  real lambda    = ((real *)params)[1];
  real phi       = ((real *)params)[2];
  real phi_match = ((real *)params)[3];
  real gamma     = ((real *)params)[4];

  real f = 0.0;

  if(x != 0.0) {
    real alpha_hat  = exp(pi*phi_match) + (1 - pow(x, gamma)) / (pow(x, gamma));
    real dalpha_hat = gamma / (pow(x, gamma+1.0)); // taken -ve to change direction of integral
    real alpha      = -exp(pi*phi);
    
    f = atan((D*sin(lambda)*pow(alpha_hat,-lambda/pi)) /
             (D*cos(lambda)*pow(alpha_hat,-lambda/pi)+1.0));
    f *= dalpha_hat;
    f /= (alpha_hat + alpha) * pi;
  }
  return f;
}

real far_upstream_integrand_s(real x, void *params)
{
  real D      = ((real *)params)[0];
  real lambda  = ((real *)params)[1];
  real phi     = ((real *)params)[2];

  real f = 0.0;

  if(x != 0.0) {
    real alpha_hat  = x;
    real dalpha_hat = 1.0;
    real alpha      = -exp(pi*phi);
    
    f = atan(-(D*sin(lambda)*pow(alpha_hat,lambda/pi)) /
              (D*cos(lambda)*pow(alpha_hat,lambda/pi)+1.0));
    f *= dalpha_hat;
    f /= (alpha_hat + alpha) * pi;
  }
  return f;
}

void removed_singularity_general_s(const vector& theta_s,
                                   const vector&     phi,
                                   const vector& phi_mid,
                                   vector&             o)
{
  unsigned int n = theta_s.length();

  for(unsigned int i = 1; i <= (n-1); i++) {
  o(i) = (theta_s(i) + theta_s(i+1)) *
          log(fabs((exp(pi*phi(n)) - exp(pi*phi_mid(i))) /
	           (exp(pi*phi(1)) - exp(pi*phi_mid(i))))) / (2.0 * pi);
  }

}

/*
real gaussian_downstream_b(real x, void *params)
{
  vector topography = vector(2);
  topography(1) = ((real *)params)[0]; // A
  real B = ((real *)params)[1];
  topography(2) = B;
  real phi = ((real *)params)[2];

  real f = 0.0;

  if(x != 0.0) {
    real gamma = B*sqrt(2); // ?
    
    real phi_hat  = (1.0-x) / (x*gamma);
    real dphi_hat = 1.0 / (gamma*x*x); // note spurious 'sign' due to integral domain direction

    real theta_b = gaussian_theta(phi_hat, topography);

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

real gaussian_upstream_b(real x, void *params)
{
  vector topography = vector(2);
  topography(1) = ((real *)params)[0]; // A
  real B = ((real *)params)[1];
  topography(2) = B;
  real phi = ((real *)params)[2];

  real f = 0.0;

  if(x != 0.0) {
    real gamma = B*sqrt(2);
    
    real phi_hat  = (x-1.0) / (x*gamma);
    real dphi_hat = 1.0 / (gamma*x*x); // integral domain direction is ok

    real theta_b = gaussian_theta(phi_hat, topography);

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
*/

real gaussian_b(real x, void *params)
{
  vector topography = vector(3);
  topography(1) = ((real *)params)[0]; // A
  real B = ((real *)params)[1];
  topography(2) = B;
  topography(3) = ((real *)params)[2]; // L
  real phi = ((real *)params)[3];

  real f = 0.0;

  if(x != 0.0) {

    real phi_hat  = atanh(x) / (B * sqrt(2) * atanh(0.5)); // ? check with mupad
    real dphi_hat = 1.0 / (B * sqrt(2) * (1.0 - x*x) * atanh(0.5)); // ? check with mupad

    real theta_b = gaussian_theta(phi_hat, topography);

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

real smoothbox_b(real x, void *params)
{
  vector topography = vector(4);

  topography(1) = ((real *)params)[0];  // A
  topography(2) = ((real *)params)[1];  // B
  topography(3) = ((real *)params)[2];  // L
  real phi_c =  ((real *)params)[3];
  topography(4) = phi_c;

  real phi   = ((real *)params)[4];

  real f = 0.0;

  if(x != 0.0) {

    real phi_hat  = phi_c * atanh(x) / atanh(0.5);
    real dphi_hat = phi_c / ((1.0 - x*x) * atanh(0.5));

    real theta_b = smoothbox_theta(phi_hat, topography);

    f = -theta_b * exp(pi*phi_hat);
    f *= dphi_hat;
    f /= exp(pi*phi_hat) + exp(pi*phi);
  }

  if(isnan(f)) f = 0;

  return f;
}

/**
  Is this game legal?
 */
void compute_convolution_symmetric_s(const int           i,
                                     const vector& theta_s,
                                     const vector&     phi,
                                     const vector&    dphi,
                                     const vector& phi_mid,
                                     vector&             o)
{
  for(int j = 1; j <= theta_s.length(); j++) {
    o(j)  = 0.5 * (theta_s(j) - ((theta_s(i) + theta_s(i+1))/2.0)) / tanh(pi*(phi(j) - phi_mid(i)) / 2.0);
    o(j) += 0.5 * (theta_s(j) + ((theta_s(i) + theta_s(i+1))/2.0)) / tanh(pi*(phi(j) + phi_mid(i)) / 2.0);
    o(j) *= dphi(j);
  }
}


real far_downstream_symmetric_integrand_s(real x, void *params)
{
  real D         = ((real *)params)[0];
  real lambda    = ((real *)params)[1];
  real phi       = ((real *)params)[2];
  real phi_match = ((real *)params)[3];
  real gamma     = ((real *)params)[4];

  real f = 0.0;

  if(x != 0.0) {
    real alpha_hat  = exp(pi*phi_match) + (1 - pow(x, gamma)) / (pow(x, gamma));
    real dalpha_hat = gamma / (pow(x, gamma+1.0));

    real alpha = -exp(pi*phi);
    // Downstream -> u =  1 + D e^{-\lambda \phi} \cos \lambda,
    //               v = -De^{-\lambda \phi} \sin \lambda.
    real f1    = atan((D*sin(lambda)*pow(alpha_hat,-lambda/pi)) / 
                      (D*cos(lambda)*pow(alpha_hat,-lambda/pi) + 1.0));
    f1        *= dalpha_hat;
    f1        /= (alpha_hat + alpha) * pi;
    f = f1;

    /*alpha   = -exp(-pi*phi);    
    real f2 = atan(Dv*sin(lambda)*pow(alpha_hat,-lambda/pi) / 
                   (Du*cos(lambda)*pow(alpha_hat,-lambda/pi) - 1.0));
    f2      *= dalpha_hat;
    f2      /= (alpha_hat + alpha) * pi;
    f = (f1 + f2) / 2.0;*/
  }
  return f;
}

real far_upstream_symmetric_integrand_s(real x, void *params)
{
  real D      = ((real *)params)[0];
  real lambda = ((real *)params)[1];
  real phi    = ((real *)params)[2];

  real f = 0.0;

  if(x != 0.0) {
    real alpha_hat  = x;
    real dalpha_hat = 1.0;
    
    real alpha = -exp(pi*phi);
    real f1    = atan(-(D*sin(lambda)*pow(alpha_hat,-lambda/pi)) / 
                       (D*cos(lambda)*pow(alpha_hat,-lambda/pi) + 1.0));
    f1        *= dalpha_hat;
    f1        /= (alpha_hat + alpha) * pi;
    f = f1;

    /*alpha   = -exp(-pi*phi);
    real f2 = atan(Dv*sin(lambda)*pow(alpha_hat,lambda/pi) /
                   (1.0-Du*cos(lambda)*pow(alpha_hat,lambda/pi)));
    f2     *= dalpha_hat;
    f2     /= (alpha_hat + alpha) * pi;
    f = (f1 + f2) / 2.0;*/
  }
  return f;
}

void removed_singularity_symmetric_s(const vector& theta_s,
                                     const vector&     phi,
                                     const vector& phi_mid,
                                     vector&             o)
{
  unsigned int n = theta_s.length();

  for(unsigned int i = 1; i <= (n-1); i++) {
  o(i) = (theta_s(i) + theta_s(i+1)) *
         log(fabs(sinh(pi*(phi(n) - phi_mid(i))/2.0) / 
                  sinh(pi*(phi(n) + phi_mid(i))/2.0))) / (2.0 * pi);
  }

}

/* NEW STUFF!

 What I need
 1. GENERAL:
    - a function that calculates the angle of the topography, theta
 2. OVER FREE SURFACE:
    * Far downstream
    * Far upstream
    - Integrand in middle part
 3. OVER BOTTOM TOPOGRAPHY:
    * function as per matlab tests that computes the integral over the bottom.
 */

void topography_compute_convolution_s(const real      phi_b,
                                      const vector& theta_s,
                                      const vector&     phi,
                                      const vector&    dphi,
                                      vector&             o)
{
  for(int j = 1; j <= theta_s.length(); j++) {
    o(j)  = theta_s(j) * exp(pi*phi(j));
    o(j) /= exp(pi*phi(j)) + exp(pi*phi_b);
    o(j) *= dphi(j);
  }
}

/**
\todo - Actually write up the maths for this.
 */
void topography_compute_symmetric_convolution_s(const real      phi_b,
                                                const vector& theta_s,
                                                const vector&     phi,
                                                const vector&    dphi,
                                                vector&             o)
{
  for(int j = 1; j <= theta_s.length(); j++) {
    o(j)  = theta_s(j) * (exp(pi*phi(j))  / (exp(pi*phi(j)) + exp(pi*phi_b)) - \
			  exp(-pi*phi(j)) / (exp(-pi*phi(j)) + exp(pi*phi_b)));
    o(j) *= dphi(j);
  }
}

/* ... Mostly checked (see Research/Box_Disturbance/test_downstream_integrand.m) */
real topography_far_upstream_s(real x, void *params)
{
  real D      = ((real *)params)[0];
  real lambda  = ((real *)params)[1];
  real phi     = ((real *)params)[2];

  real f = 0.0;
  if(x != 0.0) {
    real alpha_hat  = x;
    real dalpha_hat = 1.0;
    real alpha      = exp(pi*phi);
    
    f = atan(-(D*sin(lambda)*pow(alpha_hat,lambda/pi)) /
              (D*cos(lambda)*pow(alpha_hat,lambda/pi)+1.0));
    f *= dalpha_hat;
    f /= (alpha_hat + alpha) * pi;
  }
  return f;
}


/*  Mostly checked (see Research/Box_Disturbance/test_downstream_integrand.m) */
real topography_far_downstream_s(real x, void *params)
{
  real D         = ((real *)params)[0];
  real lambda    = ((real *)params)[1];
  real phi       = ((real *)params)[2];
  real phi_match = ((real *)params)[3];
  real gamma     = ((real *)params)[4];

  real f = 0.0;

  if(x != 0.0) {
    real alpha_hat  = exp(pi*phi_match) + (1 - pow(x, gamma)) / (pow(x, gamma)); //< \todo sometimes this will lose precision?
    real dalpha_hat = gamma / (pow(x, gamma+1.0)); // taken -ve to change direction of integral
    real alpha      = exp(pi*phi);
    
    f = atan((D*sin(lambda)*pow(alpha_hat,-lambda/pi)) /
             (D*cos(lambda)*pow(alpha_hat,-lambda/pi)+1.0));
    f *= dalpha_hat;
    f /= (alpha_hat + alpha) * pi;
  }
  return f;
}

void topography_gaussian_removed_singularity_b(const real           phi,  
                                               const real         phi_m,
                                               const vector& topography,
                                               real&                  o)
{

  real theta_b = gaussian_theta(phi, topography);

  o = -theta_b * log(fabs((exp( pi*phi_m) - exp(pi*phi)) / \
                          (exp(-pi*phi_m) - exp(pi*phi)))) / pi;
}

void topography_smoothbox_removed_singularity_b(const real           phi,  
                                                const real         phi_m,
                                                const vector& topography,
                                                real&                  o)
{

  real theta_b = smoothbox_theta(phi, topography);

  o = -theta_b * log(fabs((exp( pi*phi_m) - exp(pi*phi)) / \
                          (exp(-pi*phi_m) - exp(pi*phi)))) / pi;
}

/*Should be close to o.k.? Follows from Research/BoxDisturbance/test_integral_bottom.m */

real topography_gaussian_integrand_b(real x, void *params)
{
  vector topography = vector(3);

  topography(1) = ((real *)params)[0];  // A
  topography(2) = ((real *)params)[1];  // B
  topography(3) = ((real *)params)[2];  // L

  real phi = ((real *)params)[3];

  real f = 0.0;

  real phi_hat  = x;

  real theta_b_hat = gaussian_theta(phi_hat,topography);
  real theta_b     = gaussian_theta(phi, topography);

  f = -(theta_b_hat - theta_b) * exp(pi*phi_hat);
  f /= exp(pi*phi_hat) - exp(pi*phi);

  return f;
}

/*Should be close to o.k.? Follows from Research/BoxDisturbance/test_integral_bottom.m */

real topography_smoothbox_integrand_b(real x, void *params)
{
  vector topography = vector(4);

  topography(1) = ((real *)params)[0];  // A
  topography(2) = ((real *)params)[1];  // B
  topography(3) = ((real *)params)[2];  // L
  topography(4) = ((real *)params)[3];  // phi_c

  real phi   = ((real *)params)[4];

  real f = 0.0;

  real phi_hat  = x;

  if(fabs(phi_hat - phi) > 1e-9) {
    real theta_b_hat = smoothbox_theta(phi_hat,topography);
    real theta_b     = smoothbox_theta(phi, topography);
    /*if(fabs(temp_phi - 5.0) < 1e-5) printf("theta! %.14f %.14f \n", theta_b, theta_b_hat);
    if(fabs(temp_phi - 5.0) < 1e-5) printf("phi! %.14f %.14f \n", phi, exp(pi*phi));
    if(fabs(temp_phi - 5.0) < 1e-5) printf("phi_hat! %.14f %.14f \n", phi_hat, exp(pi*phi_hat));*/

    f = -(theta_b_hat - theta_b) * exp(pi*phi_hat);
    f /= exp(pi*phi_hat) - exp(pi*phi);
  } else {
    real mini_eps = 1e-5;
    real f1, f2;

    real theta_b_hat1 = smoothbox_theta(phi_hat+mini_eps,topography);
    real theta_b_hat2 = smoothbox_theta(phi_hat-mini_eps,topography);
    real theta_b     = smoothbox_theta(phi, topography);

    f1 = -(theta_b_hat1 - theta_b) * exp(pi*(phi_hat+mini_eps));
    f1 /= exp(pi*(phi_hat+mini_eps)) - exp(pi*phi);

    f2 = -(theta_b_hat2 - theta_b) * exp(pi*(phi_hat-mini_eps));
    f2 /= exp(pi*(phi_hat-mini_eps)) - exp(pi*phi);
    f = 0.5 * (f1 + f2);
  }

  /*if(fabs(temp_phi - 5.0) < 1e-5) printf("hello! %.14f \n", f);*/

  return f;
}
