/** \file real_vector.cpp

  Contains friend functions of the vector and matrix classes
*/

#include "real_vector.h"
#include <cmath>
#include <memory>

#ifdef __VECTOR_USE_GSL
  #include <gsl/gsl_linalg.h>
#endif

#define veclua_max(a,b) a > b ? a : b
#define veclua_min(a,b) a < b ? a : b

real empty_vec = 0.0;

int solve_lin_mat(const real_matrix& A, real_vector& x, const real_vector& b)
{
    const unsigned int m = A.rows();
    const unsigned int n = A.cols();

#ifdef __VECTOR_USE_GSL

    const unsigned int k = m < n ? n : m;

    const std::unique_ptr<real[]> heap_m(new real[(n*n)+(k*n)]);
    const std::unique_ptr<real[]> heap_v(new real[n+(2*k)]);

    real_matrix U(k, n, &(heap_m[0]));
    real_matrix V(n, n, &(heap_m[k*n]));
    real_vector s(n, &(heap_v[0]));
    real_vector w(k, &(heap_v[n]));
    real_vector b0(k, &(heap_v[n+k]));

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

#else

    /* Heap for matrices needed
     * name    size
     * R       (m,n)
     * Q       (m,m) */
    const std::unique_ptr<real[]> heap_m(new real[(m*n)+(m*m)]);
    /* Heap for vectors needed
     * name   size
     * y      (m) */
    const std::unique_ptr<real[]> heap_v(new real[m]);

    real_vector y(m,&(heap_v[0]));

    real_matrix Q(m,m,&(heap_m[0]));
    real_matrix R(m,n,&(heap_m[m*m]));

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

#endif

    return 0;

}

/** Perform QR factorisation
  
  \param A \f$m\times n\f$ input matrix
  \param Q \f$m \times m\f$ orthogonal output matrix
  \param R \f$m \times n\f$ upper-right triangular output matrix

  Uses Givens rotations to calculate the QR factorisation of a matrix.
*/
int qr_givens(const real_matrix &A, real_matrix &Q, real_matrix &R)
{
    const unsigned int m = A.rows();
    const unsigned int n = A.cols();

#ifdef __VECTOR_USE_GSL
    const unsigned int k = m < n ? m : n;

    const std::unique_ptr<real[]> heap_v(new real[k]);
    const std::unique_ptr<real[]> heap_m(new real[m*n]);

    real_matrix QR(m, n, &(heap_m[0]));
    real_vector tau(k, &(heap_v[0]));

    for(unsigned int i = 1; i <= m; i++) {
        for(unsigned int j = 1; j <= n; j++) QR(i,j) = A(i,j);
    }

    gsl_linalg_QR_decomp(QR.gsl(), tau.gsl());
    gsl_linalg_QR_unpack(QR.gsl(), tau.gsl(), Q.gsl(), R.gsl());

#else

    for(unsigned int i = 1; i <= m; i++) {
        for(unsigned int j = 1; j <= n; j++) R(i,j) = A(i,j);
        for(unsigned int j = 1; j <= m; j++) Q(i,j) = 0.0;
    }

    for(unsigned int i = 1; i <= m; i++) Q(i,i) = 1.0;

    for(unsigned int i = 1; i <= m; i++) {

        for(unsigned int k = i+1; k <= m; k++) {

            const real a = R(i,i);
            const real b = R(k,i);

            if(b != 0.0) {
                real c, s;

                if(a == 0.0) {
                    c = 0.0;
                    s = rsign(b);
                } else if(fabs(b) > fabs(a)) {
                    const real t = a / b;
                    const real u = rsign(b) * sqrt(1.0+(t*t));
                    s = 1.0 / u;
                    c = s * t;
                } else {
                    const real t = b / a;
                    const real u = rsign(a) * sqrt(1.0+(t*t));
                    c = 1.0 / u;
                    s = c * t;
                }

                for(unsigned int j = 1; j <= n; j++) {
                    const real r1 = R(i,j);
                    const real r2 = R(k,j);
                    R(i,j) = (c * r1) + (s * r2);
                    R(k,j) = (c * r2) - (s * r1);
                }

                for(unsigned int j = 1; j <= m; j++) {
                    const real q1 = Q(i,j);
                    const real q2 = Q(k,j);
                    Q(i,j) = (c * q1) + (s * q2);
                    Q(k,j) = (c * q2) - (s * q1);
                }

            } 
        }
    }

    // Q = Q'
    for(unsigned int i = 1; i <= m; i++) {
        for(unsigned int j = i+1; j <= m; j++) {
            const real t = Q(i,j);
            Q(i,j) = Q(j,i);
            Q(j,i) = t;
        }
    }

#endif

}

void veclua_pushtable(lua_State* L, const real_vector& v) {

    const int NARG = lua_gettop(L);

    lua_checkstack(L, 3);

    lua_newtable(L);
  
    for(unsigned int i = 1; i <= v.length(); i++) {
        lua_pushnumber(L, i);
        lua_pushnumber(L, v(i));
        lua_rawset(L, -3);
    }

    if(lua_gettop(L) != NARG+1) {
        lua_pushstring(L, "pushtable() failed stack-size exit condition.");
        lua_error(L);
    }

    return;

}

int veclua_veclength(lua_State* L,
                     const int index) {

    const int NARG = lua_gettop(L);

    lua_checkstack(L, 2);

    int m = 0;

    switch(lua_type(L, index)) {
        case LUA_TNUMBER:
            if(lua_gettop(L) != NARG) {
                lua_pushstring(L,
                               "veclua_veclength() failed stack-size exit "
                               "condition.");
                lua_error(L);
            }
            return 1;
        case LUA_TTABLE:
            lua_pushnil(L);
            while(lua_next(L, (index > 0) ? index : index-1)) {
                if(lua_type(L, index) != LUA_TNUMBER &&
                        lua_tonumber(L, 1) == m + 1) {
                    lua_pop(L, 1);
                    lua_pushstring(L,
                                   "veclua_veclength() requires vector stored "
                                   "with consecutive integer indices from 1 "
                                   "to n.");
                    lua_error(L);
                }
                m++;
                lua_pop(L, 1);
            }
            if(lua_gettop(L) != NARG) {
                lua_pushstring(L,
                               "veclua_veclength() failed stack-size exit "
                               "condition.");
                lua_error(L);
            }
            return m;
        default:
            return 0;
    }
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
void veclua_tovector(lua_State*    L, /**< The Lua context/stack */
                     const int index, /**< The index into Lua stack */
                     real_vector&  v  /**< Vector to modify */) {
  
    const int NARG = lua_gettop(L);

    lua_checkstack(L, 2);

    int i = 0;
  
    switch(lua_type(L, index)) {
        case LUA_TNUMBER:
            v(1) = lua_tonumber(L, index);
            return;
        case LUA_TTABLE:
            lua_pushnil(L);
            while(lua_next(L, (index > 0) ? index : index-1)) {
                v(++i) = (real)lua_tonumber(L, -1);
                lua_pop(L, 1);
            }
            if(lua_gettop(L) != NARG) {
                lua_pushstring(L,
                               "veclua_tovector() failed stack-size exit "
                               "condition.");
                lua_error(L);
            }
            return;
        default:
            return;
    }

}
 
