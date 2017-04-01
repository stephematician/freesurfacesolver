#include "integrals_pressure_formulation.h"


void compute_convolution_s(const int           i,
                           const vector&     v_s,
                           const vector&     phi,
                           const vector&    dphi,
                           const vector& phi_mid,
                           vector&             o)
{
  for(int j = 1; j <= v_s.length(); j++) {
    o(j) = (v_s(j) - ((v_s(i) + v_s(i+1))/2.0)) * exp(pi*phi(j));
    o(j) /= exp(pi*phi(j)) - exp(pi*phi_mid(i));
    o(j) *= dphi(j);
  }
}

/* phi_match should be negative */
void compute_far_upstream_integral_s(const vector& phi_mid,
                                     const real  phi_match,
                                     const real         Dv,
                                     const real     lambda,
                                     vector&             o)
{
  unsigned int n = phi_mid.length();

  for(unsigned int i = 1; i <= n; i++) {
    o(i) = -Dv * sin(lambda) * exp(-pi*(phi_mid(i) - phi_match)) * exp(lambda*phi_match) * \
            gsl_sf_hyperg_2F1(1.0, 
                              1.0+lambda/pi,
                              2.0+lambda/pi,
                              exp(-pi*(phi_mid(i) - phi_match))) / \
            (pi + lambda);
  }

}

/* phi_match should be positive */
void compute_far_downstream_integral_s(const vector& phi_mid,
                                       const real  phi_match,
                                       const real         Dv,
                                       const real     lambda,
                                       vector&             o)
{
  unsigned int n = phi_mid.length();

  for(unsigned int i = 1; i <= n; i++) {
    o(i) = -Dv * sin(lambda) * exp(-lambda*phi_match) * \
            gsl_sf_hyperg_2F1(1.0,
                              lambda/pi,
                              1.0+lambda/pi,
                              exp(-pi*(phi_match - phi_mid(i)))) / \
            lambda;
  }

}

void removed_singularity_general_s(const vector&     v_s,
                                   const vector&     phi,
                                   const vector& phi_mid,
                                   vector&             o)
{
  unsigned int n = v_s.length();

  for(unsigned int i = 1; i <= (n-1); i++) {
  o(i) = (v_s(i) + v_s(i+1)) *
          log(fabs((exp(pi*phi(n)) - exp(pi*phi_mid(i))) /
	           (exp(pi*phi(1)) - exp(pi*phi_mid(i))))) / (2.0 * pi);
  }
}

void compute_convolution_symmetric_s(const int           i,
                                     const vector&     v_s,
                                     const vector&     phi,
                                     const vector&    dphi,
                                     const vector& phi_mid,
                                     vector&             o)
{

  for(int j = 1; j <= v_s.length(); j++) {
    o(j)  = 0.5 * (v_s(j) + ((v_s(i) + v_s(i+1))/2.0)) / tanh(pi*(phi(j) + phi_mid(i))/2.0);
    o(j) += 0.5 * (v_s(j) - ((v_s(i) + v_s(i+1))/2.0)) / tanh(pi*(phi(j) - phi_mid(i))/2.0);
    o(j) *= dphi(j);
  }

}

void removed_singularity_symmetric_s(const vector&     v_s,
                                     const vector&     phi,
                                     const vector& phi_mid,
                                     vector&             o)
{
  unsigned int n = v_s.length();

  for(unsigned int i = 1; i <= (n-1); i++) {
    o(i) = (v_s(i) + v_s(i+1)) *
         log(fabs(sinh(pi*(phi(n) - phi_mid(i))/2.0) / 
                  sinh(pi*(phi(n) + phi_mid(i))/2.0))) / (2*pi);
  }

}
