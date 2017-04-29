/** \file vector.h

  A mathematical library for vector and matrix computations which interfaces
  with Lua. Provides a simple linear system solver and QR factorisation method.
  
  \todo Improve the Lua support for matrices and systems of equations.
*/

#ifndef __VECTOR_H_
#define __VECTOR_H_

#ifdef __VECTOR_USE_GSL
  #include <gsl/gsl_math.h>
  #include <gsl/gsl_matrix.h>
  #include <gsl/gsl_sf_hyperg.h>
#endif

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}

#include <cmath>

#define rsign(a) ((real)a > (real)0.0 ? 1.0 : ((real)a < (real)0.0 ? -1.0 : 0.0)) /**< Similar to copysign(1, a) except returns 0 if argument is 0. */

class real_vector;

// Must use double for GSL
#ifdef __VECTOR_USE_GSL
typedef double real;
#else
typedef double real;
#endif

extern real empty_vec;

/* These are for constructing vectors from top-of-the-stack tables */
extern void veclua_pushtable(lua_State*, const real_vector&);
extern void  veclua_tovector(lua_State*, const int, real_vector&);
extern int  veclua_veclength(lua_State*, const int);

/** A real number vector which is based on shallow copy of the data, it should
    never 'own' the data.

    Typical usage would be to allocate some space on the heap, then call the
    constructor vector(n, address-to-heap)

    Functions for pushing to or from the Lua stack utilises deep copying.
  
    Depending upon compiler option, this supports the use of GSL by providing a
    GSL instance of the vector.
*/
class real_vector {
    public:
        /** Default constructor
            
            If using GSL sets up a GSL instance with no block or owner.
        */
        real_vector(
            const unsigned int in_n = 0,      /**< number of elements in vector data */
            real*           in_data = nullptr, /**< pointer to vector data */
            lua_State*            L = nullptr, /**< lua instance containing vector */
            const int         index = 0        /**< index of vector on stack */
        ) : data((in_n == 0) ? &empty_vec : in_data),
            n_data(in_n) {

                if(L != nullptr && index != 0) {

                      const int NARG = lua_gettop(L);
                      lua_checkstack(L, 2);
                      int i = 0;
  
                      switch(lua_type(L, index)) {
                          case LUA_TNUMBER:
                              in_data[0] = lua_tonumber(L, index);
                          case LUA_TTABLE:
                              lua_pushnil(L);
                              while(lua_next(L, (index > 0) ?
                                        index :
                                        index-1)) {
                                  if(lua_type(L, -1) != LUA_TNUMBER) {
                                      lua_pop(L, 2); 
                                      lua_pushstring(
                                          L,
                                          "vector() not compatible with item "
                                          "on stack."
                                      );
                                      lua_error(L);
                                      break;
                                  }
                                  if(in_n <= i) {
                                      lua_pop(L, 2);
                                      lua_pushstring(
                                          L,
                                          "vector() not large enough for item "
                                          "on stack"
                                      );
                                      lua_error(L);
                                      break;
                                  }
                                  in_data[i++] = (real)lua_tonumber(L, -1);
                                  lua_pop(L, 1);
                              }
                          default:
                              lua_pushstring(
                                  L,
                                  "vector() not compatible with item on stack."
                              );
                              lua_error(L);
                      }
                          if(lua_gettop(L) != NARG) {
                              lua_pushstring(
                                  L,
                                  "vector() failed stack-size exit condition."
                              );
                              lua_error(L);
                          }
                }

#ifdef __VECTOR_USE_GSL
            gsl_inst.size = (in_n == 0) ? 1 : in_n;
            gsl_inst.stride = 1;
            gsl_inst.data = (in_n == 0) ? &empty_vec : data;
            gsl_inst.owner = 0;
            gsl_inst.block = 0;
#endif

        }
 
        /** Access \f$n\f$th element of the vector */
        real  operator()(const int n) const { return data[n-1]; }
        
        /** Access \f$n\f$th element of the vector */
        real& operator()(const int n) { return data[n-1]; }

#ifdef __VECTOR_USE_GSL
        /** Return a GSL instance */
        gsl_vector* gsl() { return &gsl_inst; }
        
        /** Return a GSL instance */
        const gsl_vector* gsl() const { return &gsl_inst; }
#endif

        /** Determine if any element of the vector is undefined or infinite */
        bool issingular(void) const
        {
            for(unsigned int i = 1; i <= length(); i++) {
                if(isnan(data[i-1]) || isinf(data[i-1])) return true;
            }

            return false;
        }

        /** Accessor for length of vector */
        unsigned int length() const { return n_data; }

    protected:
        real           *data; /**< pointer to vector data */
 
    private:
        unsigned int  n_data; /**< length of vector */
        bool      self_alloc;
#ifdef __VECTOR_USE_GSL
        gsl_vector gsl_inst; /**< GSL instance of the vector */
#endif

};


/** Mathematical matrix class

    A real number matrix which is based on shallow copy of the data, it should
    never 'own' the data.

    Typical usage would be to allocate some space on the heap, then call the
    constructor real_matrix(n, address-to-heap)


    

    \todo Unsure of how to support Lua for this class
*/
class real_matrix {
    public:
        /** Default constructor */
        real_matrix(const unsigned int in_m = 1,
               const unsigned int in_n = 1,
               real*           in_data = 0) :
               m_dim(in_data == 0 ? 1 : in_m),
               n_dim(in_data == 0 ? 1 : in_n),
               data(in_data == 0 ? &empty_vec : in_data) {

#ifdef __VECTOR_USE_GSL
            gsl_inst.size1 = in_data == 0 ? 1 : in_m;
            gsl_inst.size2 = in_data == 0 ? 1 : in_n;
            gsl_inst.tda = gsl_inst.size2;
            gsl_inst.data = in_data == 0 ? &empty_vec : in_data;
            gsl_inst.owner = 0;
            gsl_inst.block = 0;
#endif

        }

        /** Destructor */ 
        ~real_matrix() { }

        /** Accessor to matrix data
            @param i Row index
            @param j Column index
        */
        real  operator()(const int i, const int j) const { 
            return data[(i-1)*n_dim + (j-1)];
        }

        /** Accessor to matrix data
            @param i Row index
            @param j Column index
        */
        real&  operator()(const int i, const int j) { 
            return data[(i-1)*n_dim + (j-1)];
        }

#ifdef __VECTOR_USE_GSL
        /** Accessor for associated GNU Scientific Library data structure */
        gsl_matrix* gsl() { return &gsl_inst; }
#endif

        /** Accessor for number of rows */
        unsigned int rows() const { return m_dim; }
        /** Accessor for number of columns */
        unsigned int cols() const { return n_dim; }

        /** Solves a system of linear equations

          \param A \f$m\times n\f$ matrix of coefficients
          \param x \f$n\f$-dimensional vector to store solution in
          \param b \f$m\f$-dimensional vector for right hand side

          Uses QR factorisation method to solve a system of equations
          \f$A\mathbf{x} = \mathbf{b}\f$ by calculating
          \f[ R\mathbf{x} = Q^T \mathbf{b} \f]
          in two steps, first by computing the right hand side, then using
          back substitution to calculate \f$\mathbf{x}\f$.
        */
        friend int solve_lin_mat(const real_matrix&,
                                 real_vector&,
                                 const real_vector&);
        friend int qr_givens(const real_matrix&,
                             real_matrix&,
                             real_matrix&);

    private:
        unsigned int m_dim; /**< number of rows */
        unsigned int n_dim; /**< number of columns */
        real *data;         /**< pointer to matrix data in row-major order */ 
#ifdef __VECTOR_USE_GSL
        gsl_matrix gsl_inst; /**< Associated GNU Scientific Library data structure */
#endif

};

#endif

