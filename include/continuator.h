/** \file continuator.h

  A header file which provides functionality for numerically continuing
  solutions to systems of equations (using a predictor/corrector method).
*/

#ifndef __CONTINUATOR_H
#define __CONTINUATOR_H

#include "vector.h"
#include <math.h>
#include <new>

using namespace std;

#define PCC_DEFAULT_DELTA_X    1E-9
#define PCC_DEFAULT_STEPSIZE   1E-1
#define PCC_DEFAULT_ERR_NEWTON 1E-9
#define PCC_MAX_ATTEMPTS       10

extern "C" int luaopen_continuator(lua_State *);

/** Enumerate potential results returned by continuation process */
typedef enum {
  PROGRESS_CONVERGED,
  PROGRESS_BIFURCATION,
  PROGRESS_DIVERGED,
  PROGRESS_ERROR
} prog_result;

/** Enumerate directions of travel for continuation */
typedef enum {
  PROG_FORWARD,
  PROG_BACKWARD
} nc_dirn;

using namespace std;

/** Pure-virtual class for a contiuation method */
class continuator {
  public:
    /** Default constructor */
    continuator(const vecfun in_H = 0) : H(in_H) { }/** Numerically continue/progress solution
      
      \param x \f$n+1\f$-dimension input vector (original solution)
      \param o \f$n\f$-dimension output vector (one dimension in x is a 
               parameter)
      \param args extra arguments to pass to the governing equation function H
      \param dir direction of continuation (forward or backward)

      Progresses a solution further along an arc determined by a free
      parameter in the system.
    */

    /** Virtual destructor */
    virtual ~continuator() { }

    /** Numerically continue/progress solution
      
      \param x \f$n+1\f$-dimension input vector (original solution)
      \param tangent \f$n+1\f$-dimension tangent vector
      \param o \f$n+1\f$-dimension output vector (one dimension in x is a 
               parameter)
      \param args extra arguments to pass to the governing equation function H
      \param dir direction of continuation (forward or backward)

      Progresses a solution further along an arc determined by a free
      parameter in the system.
    */
    virtual prog_result progress(const vector& x,
                                 vector&       o,
                                 vector& tangent,
                                 void*      args,
                                 nc_dirn&    dir) const = 0;
    
    
    virtual prog_result newton(const vector& x,
                               vector&       o,
                               void*      args) const = 0;

  protected:
    vecfun H;  /**< governing set of equations */
};

/** Predictor-Corrector continuation method (E L Allgower & K Georg) */
class PCcontinuator : public continuator {
  public:
    /** Default constructor */
    PCcontinuator(const vecfun in_H = 0,
                  const real stepsize = PCC_DEFAULT_STEPSIZE)
                 : continuator(in_H),
                   h(stepsize)
    { }

    /** Destructor method */
    virtual ~PCcontinuator() { }

    /** Numerically continue/progress solution
      
      \param x \f$n+1\f$-dimension input vector (original solution)
      \param tangent \f$n+1\f$-dimension tangent vector.
      \param o \f$n+1\f$-dimension output vector
      \param args extra arguments to pass to the governing equation function H
      \param dir direction of continuation (forward or backward)

      Progresses a solution further along an arc determined by a free
      parameter in the system, using the direction induced by the
      Jacobian as a 'predictor' then a Newton method to 'correct' the
      solution. Tangent vector is used to determine the current direction of the
      system, and for detecting bifurcations.
    */
    virtual prog_result progress(const vector&,
                                 vector&,
                                 vector&,
                                 void*,
                                 nc_dirn&)
    const;
    
    /** Find solution holding continuation variable fixed using Newton's method.
      
      \param x \f$n+1\f$-dimension input vector (original solution)
      \param o \f$n+1\f$-dimension output vector
      \param args extra arguments to pass to the governing equation function H

      Fairly straightforward, holds the continuation variable constant by
      introducing a new equation to the system. Then uses QR factorisation and
      Newton's method to find solutions iteratively.
    */
    virtual prog_result newton(const vector&,
                               vector&,
                               void*)
    const;

    /** Accessor function for stepsize of predictor */
    real get_stepsize() const { return h; }
    /** Accessor function for stepsize of predictor */
    void set_stepsize(const real stepsize)
    {
      h = stepsize;
    }

  protected:
    real h; /**< step size parameter for predictor */

};

#endif
