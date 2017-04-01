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

#include <math.h>
#include <exception>

// #include <string.h>
// #include <iostream>
//
// using namespace std;

#define rsign(a) ((real)a > (real)0.0 ? 1.0 : ((real)a < (real)0.0 ? -1.0 : 0.0)) /**< Similar to copysign(1, a) except returns 0 if argument is 0. */

class vector;

typedef void (*vecfun)(const vector&, vector&, void*);

// Must use double for GSL
#ifdef __VECTOR_USE_GSL
typedef double real;
#else
typedef double real;
#endif

extern real empty_vec;

/* These are for constructing vectors from top-of-the-stack tables */
extern void  veclua_pushtable(lua_State*, const vector&);
extern vector veclua_tovector(lua_State*, const int);

extern void veclua_next(lua_State*, const int, const int);
extern void veclua_rawset(lua_State*, const int, const int);

/** A class which uses a 'soft' copy of data sent to it. In other words
  normally a heap is allocated to store the data for a vector, and
  a pointer into that heap is given to the constructor of
  the vector class, as well as the length of the vector.
  
  The vector may be operated on directly via Lua by simply being pushed onto
  the stack as userdata with appropriate accessor functions. Similarly a table
  on the top of the stack in Lua may be converted to a vector in C context by
  a friendly method.
  
  Depending upon compiler option, this supports the use of GSL by providing a
  GSL instance of the vector.
*/
class vector {
  public:
  /** Default constructor
  
  If vector data is provided, then it uses a soft-copy of the data, otherwise it
  will allocate a new array of data.
  */
  vector(const unsigned int in_n = 0, /**< number of elements in vector data */
         real*           in_data = 0  /**< pointer to vector data */) :
           data((in_n == 0) ? &empty_vec
                            : ((in_data == 0) ? new real[in_n] : in_data)),
           n_data(in_n),
           self_alloc((in_n == 0) ? false
                                  : ((in_data == 0) ? true : false))
  {
#ifdef __VECTOR_USE_GSL
    gsl_inst.size = (in_n == 0) ? 1 : in_n;
    gsl_inst.stride = 1;
    gsl_inst.data = (in_n == 0) ? &empty_vec : data;
    gsl_inst.owner = 0;
    gsl_inst.block = 0;
#endif
  }
  
  /** Copy constructor
  
  \todo fix this in the case that the original vector is empty
  
  A hard copy of the data in the original vector */
  vector(const vector& cv) :
           data(cv.n_data == 0 ? &empty_vec : new real[cv.n_data]),
           n_data(cv.n_data),
           self_alloc(cv.n_data == 0 ? false : true)
  {
    for(unsigned int i = 1; i <= n_data; i++)
      (*this)(i) = cv(i);
  }
  
  /* deep copy */
  vector& operator=(const vector& cv)
  {
    if(&cv != this) {
      if((n_data == 0) && (cv.n_data != 0)) {
        self_alloc = true;
        data = new real[cv.n_data];
        n_data = cv.n_data;
      } else {
        if(self_alloc) {
          delete [] data;
          data = new real[cv.n_data];
          n_data = cv.n_data;
        } else {
          if(n_data != cv.n_data) {
          /** \todo Throw an exception properly? */
            char errstrr[100] = "Incompatible vector assignment\n";
            throw errstrr;
          }
        }
      }
      for(unsigned int i = 0; i < n_data; i++) data[i] = cv.data[i];
    }
    return *this;
  }

  /** Destructor */
  ~vector() {
    if(self_alloc) delete [] data;
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

  A class which uses a 'soft' copy of data sent to it. In other words
  normally a heap is allocated to store the data for a matrix in row
  major order (ie each row is stored in a continuous block), and
  a pointer into that heap is given to the constructor of
  the matrix class, as well as the dimensions of the matrix.
  
  \todo not sure how to support this class in lua.
*/
class matrix {
  public:
  /** Default constructor */
  matrix(const unsigned int in_m = 1,
         const unsigned int in_n = 1,
         real*           in_data = 0) :
           m_dim(in_data == 0 ? 1 : in_m),
           n_dim(in_data == 0 ? 1 : in_n),
           data(in_data == 0 ? &empty_vec : in_data) 
  {
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
  ~matrix() { }

  /** Accessor to matrix data
    @param i Row index
    @param j Column index
  */
  real  operator()(const int i, const int j) const 
  { 
    return data[(i-1)*n_dim + (j-1)];
  }

  /** Accessor to matrix data
    @param i Row index
    @param j Column index
  */
  real&  operator()(const int i, const int j)
  { 
    return data[(i-1)*n_dim + (j-1)];
  }

#ifdef __VECTOR_USE_GSL
  /** Accessor for associated GNU Scientific Library data structure */
  gsl_matrix* gsl() { return &gsl_inst; }
#endif

/*  void print() const {
      for(unsigned int i = 0; i < m_dim; i++) {
        for(unsigned int j = 0; j < n_dim; j++) cout << data[(i*n_dim)+j] << " ";
      cout << endl;
    }
  }*/

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
  friend int solve_lin_mat(const matrix&,
                           vector&,
                           const vector&);
  friend int qr_givens(const matrix&,
                       matrix&,
                       matrix&);
  private:
  unsigned int m_dim; /**< number of rows */
  unsigned int n_dim; /**< number of columns */
  real *data;         /**< pointer to matrix data in row-major order */ 
#ifdef __VECTOR_USE_GSL
  gsl_matrix gsl_inst; /**< Associated GNU Scientific Library data structure */
#endif
};

#endif
