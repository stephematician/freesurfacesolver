#include "integrals_pressure_formulation.h"


void compute_convolution_s(const int                 i,
                           const real_vector &     v_s,
                           const real_vector &     phi,
                           const real_vector &    dphi,
                           const real_vector & phi_mid,
                           real_vector &             o) {

    for(int j = 1; j <= v_s.length(); j++) {
        o(j) = (v_s(j) - ((v_s(i) + v_s(i+1))/2.0)) * exp(pi*phi(j));
        o(j) /= exp(pi*phi(j)) - exp(pi*phi_mid(i));
        o(j) *= dphi(j);
    }

}


void compute_convolution_symmetric_s(const int                 i,
                                     const real_vector &     v_s,
                                     const real_vector &     phi,
                                     const real_vector &    dphi,
                                     const real_vector & phi_mid,
                                     real_vector &             o) {

    for(int j = 1; j <= v_s.length(); j++) {
        o(j)  = 0.5 * (v_s(j) + ((v_s(i) + v_s(i+1)) / 2.0)) /
                    tanh(pi*(phi(j) + phi_mid(i)) / 2.0);
        o(j) += 0.5 * (v_s(j) - ((v_s(i) + v_s(i+1)) / 2.0)) /
                    tanh(pi*(phi(j) - phi_mid(i)) / 2.0);
        o(j) *= dphi(j);
    }

}


/* phi_match should be negative */
void far_upstream_integral_s(const real_vector & phi_mid,
                             const real        phi_match,
                             const real               Dv,
                             const real           lambda,
                             real_vector &             o) {

    const unsigned int n = phi_mid.length();

    for(unsigned int i = 1; i <= n; i++) {
        o(i) = -Dv * sin(lambda) * 
                    exp(-pi*(phi_mid(i) - phi_match)) * exp(lambda*phi_match) * 
                    gsl_sf_hyperg_2F1(
                        1.0, 
                        1.0+lambda/pi,
                        2.0+lambda/pi,
                        exp(-pi*(phi_mid(i) - phi_match))
                    ) / (pi + lambda);
    }

}


/* phi_match should be positive */
void far_downstream_integral_s(const real_vector & phi_mid,
                               const real        phi_match,
                               const real               Dv,
                               const real           lambda,
                               real_vector&              o) {

    const unsigned int n = phi_mid.length();

    for(unsigned int i = 1; i <= n; i++) {
        o(i) = -Dv * sin(lambda) * exp(-lambda*phi_match) *
                    gsl_sf_hyperg_2F1(
                        1.0,
                        lambda/pi,
                        1.0+lambda/pi,
                        exp(-pi*(phi_match - phi_mid(i)))
                    ) / lambda;
    }

}


void removed_singularity_s(const real_vector &     v_s,
                           const real_vector &     phi,
                           const real_vector & phi_mid,
                           real_vector &             o) {

    const unsigned int n = v_s.length();

    for(unsigned int i = 1; i <= (n-1); i++) {
        o(i) = (v_s(i) + v_s(i+1)) *
                   log(
                       fabs(
                           (exp(pi*phi(n)) - exp(pi*phi_mid(i))) /
                               (exp(pi*phi(1)) - exp(pi*phi_mid(i)))
                       )
                   ) / (2.0 * pi);
    }

}


void removed_singularity_symmetric_s(const real_vector &     v_s,
                                     const real_vector &     phi,
                                     const real_vector & phi_mid,
                                     real_vector &             o) {
    const unsigned int n = v_s.length();

    for(unsigned int i = 1; i <= (n-1); i++) {
        o(i) = (v_s(i) + v_s(i+1)) *
                   log(
                       fabs(
                           sinh(pi*(phi(n) - phi_mid(i)) / 2.0) / 
                               sinh(pi*(phi(n) + phi_mid(i)) / 2.0)
                       )
                   ) / (2.0 * pi);
    }

}

