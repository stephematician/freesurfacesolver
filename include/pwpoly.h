/** \file pwpoly.h

    Function declarations (external) for managing piecewise polynomial
    transformations of domain

*/

#ifndef __PWPOLY_H_
#define __PWPOLY_H_

#include "real_vector.h"

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}


/** \internal Transform (arbitrary) relating \f$\beta\f$ to \f$\phi\f$. */
typedef void (*transform_func)(const real_vector &, const real_vector &,
                               const real_vector &, const real_vector &,
                                     real_vector &, const real);



/** Calculates a piece-wise polynomial relationship between \f$\beta\f$ and 
    \f$\phi\f$ based on a homotopy between a linear and a cubic interpolating
    spline.

    Calculates a piece-wise polynomial relationship between two variables,
    \f$\beta\f$ and \f$\phi\f$, which is based on a homotopy parameter \f$s\f$
    which varies the relationship between a piecewise linear spline and a
    piecewise cubic spline. The genereal piecewise polynomial/homotopy is
    defined on the unit intervals of \f$\beta\f$ given by \f$[i-1,i]\f$ for
    \f$i = 1,2,\ldots,n\f$. Then the piecewise polynomial/homotopy is

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

    \param beta A \f$m\f$-dimensional vector specifying the points at which the
                polynomial will be evaluated. 
    \param beta_sub A \f$n+1\f$ dimensional vector specifying data to construct
                    piecewise polynomial.
    \param phi_sub A \f$n+1\f$-dimensional vector specifying data to construct
                   the piecewise polynomial.
    \param dphi_sub A \f$n+1\f$-dimensional vector specifying data to construct
                    the piecewise polynomial.
    \param o A \f$m\f$-dimensional vector with the output of the function
             evaluated at the input \f$\beta\f$
    \param s A homotopy parameter that should be between 0 and 1.

*/
extern void pwlinear_cubic_homotopy_phi(const real_vector &     beta,
                                        const real_vector & beta_sub,
                                        const real_vector &  phi_sub,
                                        const real_vector & dphi_sub,
                                        real_vector &              o,
                                        const real                 s);

/** Calculates a piece-wise polynomial relationship between \f$\beta\f$ and 
    \f$\phi\f$ based on a homotopy between a linear and a cubic interpolating
    spline.

    Calculates a the derivative of a piece-wise polynomial relationship between
    two variables, \f$\beta\f$ and \f$\phi\f$, which is based on a homotopy
    parameter \f$s\f$ which varies the relationship between a piecewise linear
    spline and a piecewise cubic spline. The genereal piecewise polynomial/homotopy is
    defined on the unit intervals of \f$\beta\f$ given by \f$[i-1,i]\f$ for
    \f$i = 1,2,\ldots,n\f$.

    \seealso pwlinear_cubic_homotopy_phi

    \param beta A \f$m\f$-dimensional vector specifying the points at which the
                polynomial will be evaluated. 
    \param beta_sub A \f$n+1\f$ dimensional vector specifying data to construct
                    piecewise polynomial.
    \param phi_sub A \f$n+1\f$-dimensional vector specifying data to construct
                   the piecewise polynomial.
    \param dphi_sub A \f$n+1\f$-dimensional vector specifying data to construct
                    the piecewise polynomial.
    \param o A \f$m\f$-dimensional vector with the output of the function
             evaluated at the input \f$\beta\f$
    \param s A homotopy parameter that should be between 0 and 1.

*/
extern void pwlinear_cubic_homotopy_dphi(const real_vector &     beta,
                                         const real_vector & beta_sub,
                                         const real_vector &  phi_sub,
                                         const real_vector & dphi_sub,
                                         real_vector &              o,
                                         const real                 s);
/** Lua interface to pwlinear_cubic_homotopy_phi 
 
    \todo Document properly */
extern int pwlinear_cubic_homotopy_phi_luaf(lua_State*);


/** Lua interface to pwlinear_cubic_homotopy_dphi 
 
    \todo Document properly */
extern int pwlinear_cubic_homotopy_dphi_luaf(lua_State*);


#endif

