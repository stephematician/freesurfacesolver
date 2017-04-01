/** \file solver.h
 
  Provides an interface for iterative methods for solving systems of equations.

  Defines an implementation of the interface which uses Newton's method to solve
  the systems of equations.
*/
#ifndef __SOLVER_H
#define __SOLVER_H

#include "vector.h"
#include <math.h>

#define SOLVER_MAX_ATTEMPTS 11
#define SOLVER_DEFAULT_TOL  1E-6
#define SOLVER_DEFAULT_H    1E-7

class solver; 
class Newtonsolver;

/** Enumerate possible results of convergence */
typedef enum {
  SOLVER_CONVERGED,
  SOLVER_NOTCONVERGED,
  SOLVER_DIVERGED,
  SOLVER_ERROR
} solve_result;

/** Pure virtual class for iteratively solving systems of equations. */
class solver {
  public:
    /** Default constructor
    
    \param in_H Function that calculates the residual of the system of equations
                given an approximate solution
    \param in_k The number of equations
     */
    solver(const vecfun       in_H = 0,
           const unsigned int in_k = 1) : 
             H(in_H),
             n_row(in_k) { }

    /** Virtual destructor */
    virtual ~solver() { }

    /** Improve solution by performing one iteration

      \param x \f$n\f$-dimensional input vector (initial approximation)
      \param o \f$n\f$-dimensional output vector (improved solution)
      \param args extra args to pass to governing equations function H

      Pure virtual method to improve the solution by one step using
      a numerical root-finding method.
    */
    virtual int improve(const vector& x,
                        vector&       o,
                        void*      args) const = 0;

    /** Find root by iteration

      \param x \f$n\f$-dimensional input vector (initial approximation)
      \param o \f$n\f$-dimensional output vector (improved solution)
      \param args extra args to pass to governing equations function H

      Pure virtual method to find an as-accurate-as-possible solution to
      the governing set of equations.
    */
    virtual solve_result converge(const vector& x,
                                   vector&       o,
                                   void*      args) const = 0;
    /** Find root by iteration

      \param x \f$n\f$-dimensional input vector (initial approximation)
      \param o \f$n\f$-dimensional output vector (improved solution)
      \param args extra args to pass to governing equations function H
      \param h_ accuracy of solution

      Pure virtual method to find a solution such that the \f$L_{sup}\f$ norm 
      of \f$H(x)\f$ is less than \c h_.
    */
    virtual solve_result converge(const vector& x,
                                  vector&       o,
                                  void*      args,
                                  const real   h_) const = 0;
    
  protected:
    vecfun       H;     /**< Function that calculates the residual of the system
                             of equations. Assume \f$ H : \mathcal{R}^{m} \to
                             \mathcal{R}^n\f$. */
    unsigned int n_row; /**< Dimensions of return argument for H */

};

/** Newton's method for root 'improving'

  Uses Newton's method to find improved solutions to systems of equations
  using the linear equation system solving methods supplied by vector.h
  and a finite difference approximation to Jacobians.
*/
class Newtonsolver : public solver
{
  public:
    /** Default constructor */
    Newtonsolver(const vecfun       in_H = 0,
                 const unsigned int in_k = 0) :
                   solver(in_H, in_k), 
                   h(SOLVER_DEFAULT_H) { }

    /** Virtual destructor */
    virtual ~Newtonsolver() { }

    int improve(const vector&, vector&, void*) const;
    solve_result converge(vector&, void*) const;
    solve_result converge(const vector&, vector&, void*) const;
    solve_result converge(const vector&, vector&, void*, const real) const;
    
    /** Accessor function for \c h */
    real get_h() const { return h; }
    /** Accessor function for \c h */
    void set_h(const real in_h) { h = in_h; }

  protected:
    real h; /**< Jacobian \f$J\f$ used for Newton's method is calculated via the 
                 finite difference approximation 
                 \f[ J = [J_{ij}] = \frac{\partial H_i}{\partial x_j} \approx
      \frac{H_i(\mathbf{x} + h \mathbf{e}_j) - H_i(\mathbf{x})}{h} \f] where
                 \f$H_i(\mathbf{x}) = H(\mathbf{x}) \cdot \mathbf{e}_i\f$ */
};

#endif
