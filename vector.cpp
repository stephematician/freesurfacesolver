/** \file vector.cpp

  Contains friend functions of the vector and matrix classes
*/

#include "vector.h"
#include <math.h>
#include <stdio.h>

#ifdef __VECTOR_USE_GSL
  #include <gsl/gsl_linalg.h>
#endif

#define veclua_max(a,b) a > b ? a : b
#define veclua_min(a,b) a < b ? a : b

real empty_vec = 0.0;

int solve_lin_mat(const matrix& A, vector& x, const vector& b)
{
  unsigned int m = A.rows();
  unsigned int n = A.cols();
#ifdef __VECTOR_USE_GSL
  unsigned int k = m < n ? n : m;

  real* heap_m = new real[(n*n)+(k*n)];
  real* heap_v = new real[n+k+k];

  matrix U  = matrix(k, n, &(heap_m[0]));
  matrix V  = matrix(n, n, &(heap_m[k*n]));
  vector s  = vector(n, &(heap_v[0]));
  vector w  = vector(k, &(heap_v[n]));
  vector b0 = vector(k, &(heap_v[n+k]));

  for(unsigned int i = 1; i <= k; i++) {
    if(i > m) {
      b0(i) = 0;
      for(unsigned int j = 1; j <= n; j++) U(i,j) = 0;
    } else {
      b0(i) = b(i);
      for(unsigned int j = 1; j <= n; j++) U(i,j) = A(i,j);
    }
  }

  gsl_linalg_SV_decomp(U.gsl(), V.gsl(), s.gsl(), w.gsl());
  gsl_linalg_SV_solve(U.gsl(), V.gsl(), s.gsl(), b0.gsl(), x.gsl());

  delete [] heap_m;
  delete [] heap_v;

#else
  /* Heap for matrices needed
     name    size
     R       (m,n)
     Q       (m,m)
  */
  real *heap_m = new real[(m*n)+(m*m)];
  /* Heap for vectors needed
     name   size
     y      (m)
  */
  real *heap_v = new real[m];

  vector y = vector(m,&(heap_v[0]));

  matrix Q = matrix(m,m,&(heap_m[0]));
  matrix R = matrix(m,n,&(heap_m[m*m]));

  qr_givens(A, Q, R);

  /* Matrix multiplication y = Q'b */
  for(unsigned int i = 1; i <= m; i++) {
    y(i) = 0.0;
    for(unsigned int j = 1; j <= m; j++) y(i) += Q(j,i) * b(j);
  }

  for(unsigned int i = m+1; i <= n; i++) x(i) = 0.0;

  /* Use back substitution Rx = y */
  for(unsigned int i = m; i >= 1; i--) {
    if(i <= n) {
      x(i) = y(i);
      for(unsigned int j = i+1; j <= n; j++) x(i) -= R(i,j) * x(j);
      x(i) /= R(i,i);
    }
  }

  delete [] heap_m;
  delete [] heap_v;
#endif

  return 0;
}

/** Perform QR factorisation
  
  \param A \f$m\times n\f$ input matrix
  \param Q \f$m \times m\f$ orthogonal output matrix
  \param R \f$m \times n\f$ upper-right triangular output matrix

  Uses Givens rotations to calculate the QR factorisation of a matrix.
*/
int qr_givens(const matrix &A, matrix &Q, matrix &R)
{
  unsigned int m = A.rows();
  unsigned int n = A.cols();

#ifdef __VECTOR_USE_GSL
  unsigned int k = m < n ? m : n;

  real* heap_v = new real[k];
  real* heap_m = new real[m*n];

  matrix QR  = matrix(m, n, &(heap_m[0]));
  vector tau = vector(k, &(heap_v[0]));

  for(unsigned int i = 1; i <= m; i++) {
    for(unsigned int j = 1; j <= n; j++) {
      QR(i,j) = A(i,j);
    }
  }

  gsl_linalg_QR_decomp(QR.gsl(), tau.gsl());
  gsl_linalg_QR_unpack(QR.gsl(), tau.gsl(), Q.gsl(), R.gsl());

  delete [] heap_v;
  delete [] heap_m;

#else
  real a, b, t, u;
  real c, s, r1, r2, q1, q2;

  for(unsigned int i = 1; i <= m; i++) {
    for(unsigned int j = 1; j <= n; j++) R(i,j) = A(i,j);
    for(unsigned int j = 1; j <= m; j++) Q(i,j) = 0.0;
  }

  for(unsigned int i = 1; i <= m; i++) Q(i,i) = 1.0;

  for(unsigned int i = 1; i <= m; i++) {
    for(unsigned int k = i+1; k <= m; k++) {
      a = R(i,i);
      b = R(k,i);
      if(b != 0.0) {
        if(a == 0.0) {
          c = 0.0;
          s = rsign(b);
        } else if(fabs(b) > fabs(a)) {
          t = a / b;
          u = rsign(b) * sqrt(1.0+(t*t));
          s = 1.0 / u;
          c = s * t;
        } else {
          t = b / a;
          u = rsign(a) * sqrt(1.0+(t*t));
          c = 1.0 / u;
          s = c * t;
        }
        for(unsigned int j = 1; j <= n; j++) {
          r1 = R(i,j);
          r2 = R(k,j);
          R(i,j) = (c * r1) + (s * r2);
          R(k,j) = (c * r2) - (s * r1);
        }
        for(unsigned int j = 1; j <= m; j++) {
          q1 = Q(i,j);
          q2 = Q(k,j);
          Q(i,j) = (c * q1) + (s * q2);
          Q(k,j) = (c * q2) - (s * q1);
        }
      } 
    }
  }

  // Q = Q'
  for(unsigned int i = 1; i <= m; i++) {
    for(unsigned int j = i+1; j <= m; j++) {
      t = Q(i,j);
      Q(i,j) = Q(j,i);
      Q(j,i) = t;
    }
  }
#endif
}

void veclua_pushtable(lua_State* L, const vector& v)
{
  int NARG = lua_gettop(L);
  
  lua_checkstack(L, lua_gettop(L)+3);

  lua_newtable(L);
  
  for(unsigned int i = 1; i <= v.length(); i++) {
    lua_pushnumber(L, i);
    lua_pushnumber(L, v(i));
    lua_rawset(L, -3);
  }
  if(lua_gettop(L) != NARG+1) {
    lua_pushstring(L, "veclua_pushtable exit");
    lua_error(L);
  }
  return;
}

/** Makes a vector object from a table on the top of the stack

Converts a variable in Lua to a vector when the Lua variable is compatible.
If the Lua variable is a single number, then an mvector is created which has
a single element which has the same value as the Lua number.

If the Lua variable is a table which is effectively an array with indices
\f$1,2,\ldots,N\f$ then an mvector with \f$N\f$ elements is created with the
same values (in the same order as iterating over the Lua table) as the Lua
table.

The Lua stack will remained unchanged, aside from possible resizing to
accomodate iteration over a Lua table.
*/
vector veclua_tovector(lua_State*    L, /**< The Lua context/stack */
                       const int index  /**< The index into Lua stack*/)
{
  int NARG = lua_gettop(L);

  lua_checkstack(L, lua_gettop(L)+2);
  
  switch(lua_type(L, index)) {
  case LUA_TNUMBER : {
    vector v(1);
    v(1) = lua_tonumber(L, index);
    if(lua_gettop(L) != NARG) {
      lua_pushstring(L, "veclua_tovector exit");
      lua_error(L);
    }
    return v;
  }
  break;
  case LUA_TTABLE : {
    lua_pushnil(L);
    int m = 0;
    while(lua_next(L, (index > 0) ? index : index-1)) {
      m++;
      lua_pop(L, 1);
    }
    vector v(m);
    lua_pushnil(L);
    while(lua_next(L, (index > 0) ? index : index-1)) {
      /* value
         key
         ...
         table
      */
      int i = (int)lua_tonumber(L, -2);
      real vi = (real)lua_tonumber(L, -1);
      v(i) = vi;
      lua_pop(L, 1);
    }
    if(lua_gettop(L) != NARG) {
      lua_pushstring(L, "veclua_tovector exit");
      lua_error(L);
    }
    return v;
  }
  break; 
  default :
    return 0;
  }
}
