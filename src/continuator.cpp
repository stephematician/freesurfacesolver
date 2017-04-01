#include "continuator.h"

int lua_continuator_new(lua_State*);
int lua_continuator__gc(lua_State*);
int lua_continuator__index(lua_State*);
int lua_continuator__newindex(lua_State*);
int lua_continuator_progress(lua_State*);
int lua_continuator_newton(lua_State*);
int lua_continuator_set_stepsize(lua_State*);
int lua_continuator_get_stepsize(lua_State*);
void lua_continuator_proglayer(const vector&, vector&, void*);

struct lua_constant {
  const char *str;
  int       value;
};

static struct lua_constant continuator_constants[] = {
  {"CONVERGED", PROGRESS_CONVERGED},
  {"BIFURCATION", PROGRESS_BIFURCATION},
  {"DIVERGED", PROGRESS_DIVERGED},
  {"ERROR", PROGRESS_ERROR},
  {"BACKWARD", PROG_BACKWARD},
  {"FORWARD", PROG_FORWARD},  
  {NULL, 0}
};

static const luaL_Reg continuatorR[] = {
  {"new", lua_continuator_new},
  {NULL, NULL}
};

// Add constants to the table on top of the stack
// stolen from luaglfw.c
static void add_constants(lua_State *L, struct lua_constant *cn)
{
  int NARG = lua_gettop(L);
  while(cn->str)
  {
    lua_pushstring(L, cn->str);
    lua_pushnumber(L, cn->value);
    lua_rawset(L, -3);
    ++cn;
  }
  if(lua_gettop(L) != NARG) {
    lua_pushstring(L, "add_constants exit");
    lua_error(L);
  }
}

int luaopen_continuator(lua_State *L) {
  int NARG = lua_gettop(L);
  luaL_setfuncs(L, continuatorR, 0);
  
  add_constants(L, continuator_constants);
  
  if(lua_gettop(L) != NARG+1) {
    lua_pushstring(L, "luaopen_continuator_cpp exit");
    lua_error(L);
  }
  return 1;
}

prog_result PCcontinuator::progress(const vector &x0,
                                    vector&        x,
                                    vector&  tangent,
                                    void*       args,
                                    nc_dirn&     dir) const
{
  unsigned int m = x0.length()-1;
  unsigned int n   = x0.length();
  bool progress_diverged = false;
  bool        bifurcated = false;
  unsigned int attempts      = 0;

  real jj;

  real delta_x = PCC_DEFAULT_DELTA_X; /* delta for Jacobian estimate */
  real err_tol = PCC_DEFAULT_ERR_NEWTON; /* tolerance for Newton's method */
  real err_val = 2 * err_tol;

  real dot_tangents = 0.0;
  real sign_detJ    = 1.0;
  real direction    = (dir == PROG_FORWARD) ? 1.0 : -1.0;

  /* bail if H is not defined */
  if(H == 0) return PROGRESS_ERROR;

  /* Heap for vectors needed
     name   size
     predict (n)
     x1      (n)
     b       (n)
     y0      (m)
     y1      (m)
  */
  real *heap_v = new real[(3*n)+(2*m)];

  /* Heap for matrices needed
     name    size

     J       (n,m)
     R       (n,m)
     Q       (n,n)
  */
  real *heap_m = new real[(2*m*n)+(n*n)];

  vector x1      = vector(n, &(heap_v[0]));
  vector predict = vector(n, &(heap_v[n]));
  vector b       = vector(n, &(heap_v[2*n]));
  vector y0      = vector(m, &(heap_v[3*n]));
  vector y1      = vector(m, &(heap_v[(3*n)+m]));

  matrix J       = matrix(n, m, &(heap_m[0]));
  matrix R       = matrix(n, m, &(heap_m[m*n]));
  matrix Q       = matrix(n, n, &(heap_m[m*n*2]));

  H(x0, y0, args);
    
  for(unsigned int i = 1; i <= n; i++) x1(i) = x0(i);
  /*real max_jac = -INFINITY;
  int max_jac_i, max_jac_j;*/

  /* Calculate the transpose of the Jacobian of the system at x0 */
  for(unsigned int j = 1; j <= n; j++) {
    x1(j) += delta_x;
    H(x1, y1, args);

    for(unsigned int i = 1; i <= m; i++) {
      J(j,i) = (y1(i) - y0(i)) / delta_x;
    /*  if (fabs(J(j,i)) > max_jac) {
        max_jac = fabs(J(j,i));
        max_jac_i = i;
        max_jac_j = j;
      }*/
    }
    x1(j) = x0(j);
  }

  //printf("max_jac = %f, at (%i,%i)\n", max_jac, max_jac_i, max_jac_j);
  qr_givens(J, Q, R);
/*  real max_r = -INFINITY;
  int max_r_i;
  real min_r = INFINITY;
  int min_r_i;

  for(unsigned int i = 1; i <= m; i++) {
    if (fabs(R(i,i)) > max_r) {
        max_r = fabs(R(i,i));
        max_r_i = i;
    }
    if (fabs(R(i,i)) < min_r) {
        min_r = fabs(R(i,i));
        min_r_i = i;
    }
  }

  printf("max_r = %f, at (%i,%i)\n", max_r, max_r_i, max_r_i);
  printf("min_r = %f, at (%i,%i)\n", min_r, min_r_i, min_r_i);

  printf("Q(:,:) ");
  for(unsigned int i = 1; i <= n; i++) {
    for(unsigned int j = 1; j <= n; j++) {
      printf("%f ", Q(i,j));
    }
    printf("\n");
  }
  */
  for(unsigned int i = 1; i <= m; i++) sign_detJ *= rsign(R(i,i));
  
  for(unsigned int i = 1; i <= n; i++) {
    dot_tangents += Q(i,n) * tangent(i);
    tangent(i) = Q(i,n) * sign_detJ;
  }
  
  dot_tangents *= sign_detJ;
  if(dot_tangents < 0.0) {
    direction *= -1.0;
    dir = (dir == PROG_FORWARD) ? PROG_BACKWARD : PROG_FORWARD;
    bifurcated = true;
  }

  for(unsigned int i = 1; i <= n; i++)
    predict(i) = x0(i) + h * direction * tangent(i);

  /*  predicted solution needs refinement via Newton's method */
  H(predict, y0, args);

  while(err_val >= err_tol) {
    if(attempts > PCC_MAX_ATTEMPTS) {
      progress_diverged = true;
      break;
    }

    for(unsigned int i = 1; i <= n; i++) x1(i) = predict(i);

    /* Calculate the transpose of the Jacobian of the system at predict */
    for(unsigned int j = 1; j <= n; j++) {
      x1(j) += delta_x;
      H(x1, y1, args);

      for(unsigned int i = 1; i <= m; i++) J(j,i) = (y1(i) - y0(i)) / delta_x;
      x1(j) = predict(j);
    }

    qr_givens(J, Q, R);

    for(unsigned int i = 1; i <= m; i++) {
      b(i) = -y0(i);
      for(unsigned int j = 1; j <= i-1; j++) b(i) -= R(j,i) * b(j);
      b(i) /= R(i,i);
    }

    b(n) = 0.0;

    for(unsigned int i = 1; i <= n; i++) {
      for(unsigned int j = 1; j <= n; j++) predict(i) += Q(i,j) * b(j);
    }

    H(predict, y0, args);
    err_val = 0.0;
    int po = 0;
    for(unsigned int i = 1; i <= m; i++) {err_val = err_val < fabs(y0(i)) ?
                                                             fabs(y0(i)) :
                                                                   err_val;
                                                                   po = i;}
    attempts++;
  }

  for(unsigned int i = 1; i <= n; i++) x(i) = predict(i);

  delete [] heap_v;
  delete [] heap_m;
  
  return progress_diverged ? PROGRESS_DIVERGED :
                             ( bifurcated ? PROGRESS_BIFURCATION :
                                            PROGRESS_CONVERGED );
}

prog_result PCcontinuator::newton(const vector &x0,
                                  vector&        x,
                                  void*       args) const
{
  unsigned int m = x0.length()-1;
  unsigned int n   = x0.length();
  bool progress_diverged = false;
  unsigned int attempts      = 0;

  real delta_x = PCC_DEFAULT_DELTA_X; /* delta for Jacobian estimate */
  real err_tol = PCC_DEFAULT_ERR_NEWTON; /* tolerance for Newton's method */
  real err_val = 2 * err_tol;

  real dot_tangents = 0.0;
  real sign_detJ    = 1.0;

  /* bail if H is not defined */
  if(H == 0) return PROGRESS_ERROR;

  /* Heap for vectors needed
     name   size
     x1      (n)
     b       (n)
     dx      (n)
     y0      (m)
     y1      (m)
  */
  real *heap_v = new real[(4*n)+(2*m)];

  /* Heap for matrices needed
     name    size

     J       (n,n)
     R       (n,n)
     Q       (n,n)
  */
  real *heap_m = new real[(3*n*n)];

  vector x1      = vector(n, &(heap_v[0]));
  vector soln    = vector(n, &(heap_v[n]));
  vector b       = vector(n, &(heap_v[2*n]));
  vector dx      = vector(n, &(heap_v[3*n]));
  vector y0      = vector(m, &(heap_v[4*n]));
  vector y1      = vector(m, &(heap_v[(4*n)+m]));

  matrix J       = matrix(n, n, &(heap_m[0]));
  matrix R       = matrix(n, n, &(heap_m[n*n]));
  matrix Q       = matrix(n, n, &(heap_m[n*n*2]));
  
  printf("Lua continuation proglayer call\n");
  H(x0, y0, args);
  printf("Lua continuation proglayer success\n");
    
  for(unsigned int i = 1; i <= n; i++) soln(i) = x0(i);
      
  while(err_val >= err_tol) {
    if(attempts > PCC_MAX_ATTEMPTS) {
      progress_diverged = true;
      break;
    }
    printf("Lua continuation newton iteration\n");
    for(unsigned int i = 1; i <= n; i++) x1(i) = soln(i);

    /* Calculate the Jacobian of the system at x0 */
    for(unsigned int j = 1; j <= n; j++) {
      x1(j) += delta_x;
      H(x1, y1, args);

      for(unsigned int i = 1; i <= m; i++) J(i,j) = (y1(i) - y0(i)) / delta_x;
      if(j != n) {
        J(n,j) = 0.0;
      } else {
        J(n,n) = 1.0;
      }
      x1(j) = soln(j);
    }

    qr_givens(J, Q, R);

    /* Calculate Q^t y */
    for(unsigned int i = 1; i <= n; i++) {
      b(i) = Q(n,i) * (soln(n) - x0(n));
      for(unsigned int j = 1; j <= m; j++) b(i) += Q(j,i) * y0(j);
    }

    /* Solve R \delta x = b */
    /*printf("dx =");*/
    dx(n) = -b(n) / R(n,n);
    /*printf(" %f", dx(n));*/
    for(unsigned int i = m; i >= 1; i--) {
      dx(i) = -b(i);
      for(unsigned int j = i+1; j <= n; j++) dx(i) -= R(i,j) * dx(j);
      dx(i) /= R(i,i);
      /*printf(" %f", dx(i));*/
    }
    /*printf("\n");*/

    for(unsigned int i = 1; i <= n; i++)  soln(i) += dx(i);
    
    printf("Lua continuation newton step complete");
    /*printf(", here's x0 and y0\n");
    printf("x0 =");
    for(unsigned int i = 1; i <= n; i++) {
      printf(" %f", soln(i));
    }
    printf("\n");*/
    
    H(soln, y0, args);
    
    /*printf("y0 =");
    for(unsigned int i = 1; i <= m; i++) {
      printf(" %f", y0(i));
    }
    printf(" %f", soln(n) - x0(n));
    printf("\n");*/
    err_val = 0.0;
    int po = 0;
    for(unsigned int i = 1; i <= m; i++) err_val = err_val < fabs(y0(i)) ?
                                                             fabs(y0(i)) :
                                                                  err_val;

    err_val = err_val < fabs(soln(n) - x0(n)) ?
                        fabs(soln(n) - x0(n)) :
                        err_val;

    attempts++;
  }

  for(unsigned int i = 1; i <= n; i++) x(i) = soln(i);

  delete [] heap_v;
  delete [] heap_m;
  
  return progress_diverged ? PROGRESS_DIVERGED : PROGRESS_CONVERGED;
}


/** Creates a lua table which provides an interface to a PCcontinuator object.

  \param L Lua state
  
  Called from Lua, this function requires one argument - which is a valid
  residual function for continuation.
  
  A residual function takes the form (in lua) of 
  \code
  residual = function(unknowns, ...)
  \endcode
  where \c residual is a table with numerical values only (i.e. a vector).

*/
int lua_continuator_new(lua_State* L)
{
  static const int NARG = 1;
  
  if(lua_gettop(L) != 1) {
    lua_pushstring(L, "lua_continuator_new entry");
    lua_error(L);
  }
  
  switch(lua_type(L, 1))
  {
  case LUA_TFUNCTION :
    lua_newtable(L);

    lua_pushstring(L, "residual");
    lua_pushvalue(L, 1);
    lua_rawset(L, -3); // set residual to first argument
    
    lua_remove(L, 1);
    
    lua_pushstring(L, "get_stepsize");
    lua_pushcfunction(L, lua_continuator_get_stepsize);
    lua_rawset(L, -3);
    
    lua_pushstring(L, "set_stepsize");
    lua_pushcfunction(L, lua_continuator_set_stepsize);
    lua_rawset(L, -3);
    
    lua_pushstring(L, "newton");
    lua_pushcfunction(L, lua_continuator_newton);
    lua_rawset(L, -3);
    
    lua_pushstring(L, "progress");
    lua_pushcfunction(L, lua_continuator_progress);
    lua_rawset(L, -3);
    break;
  default:
    lua_pop(L, NARG);
    return 0;
  }

  lua_pushstring(L, "ud");  // make PCcontinuatour w/- custom metatable

  void* buffer = lua_newuserdata(L, sizeof(PCcontinuator));
  
  PCcontinuator* nc = new(buffer) PCcontinuator(lua_continuator_proglayer);
  

  lua_newtable(L);   // Make new metatable for PCcontinuator
  
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, lua_continuator__gc);
  lua_rawset(L, -3);

  lua_pushstring(L, "__newindex");
  lua_pushcfunction(L, lua_continuator__newindex);
  lua_rawset(L, -3);

  lua_pushstring(L, "__index");
  lua_pushcfunction(L, lua_continuator__index);
  lua_rawset(L, -3);

  lua_pushstring(L, "__metatable");
  lua_pushvalue(L, -2);
  lua_rawset(L, -3);
  
  lua_setmetatable(L, -2); // sets PCcontinuator's metadata
  
  lua_rawset(L, -3); // sets 'ud' to newuserdata
  
  if(lua_gettop(L) != 1) {
    lua_pushstring(L, "lua_continuator_new exit");
    lua_error(L);
  }
  
  return 1;
}

void lua_continuator_proglayer(const vector& input, vector& output, void* args)
{
  lua_State* L = (lua_State*)args;
  
  const int NARG = lua_gettop(L);
  
  /* Basically take the function and its arguments, and copy them onto the
  top of the stack, so as to call the function yet preserve the stack. */
  lua_pushvalue(L, 1);
  veclua_pushtable(L, input);
  for(unsigned int i = 2; i <= NARG; i++) {
    lua_pushvalue(L, i);
  }
  lua_call(L, NARG, 1);

  output = veclua_tovector(L, -1);
  lua_pop(L, 1);
  
  if(NARG != lua_gettop(L)) {
    lua_pushstring(L, "lua_continuator_proglayer exit");
    lua_error(L);
  }
}

int lua_continuator_set_stepsize(lua_State* L)
{
  if(lua_gettop(L) != 2) {
    lua_pushstring(L, "lua_continuator_set_stepsize entry");
    lua_error(L);
  }
  real _ssize = lua_tonumber(L, 2);
  lua_pop(L, 1);
  lua_pushstring(L, "ud");
  lua_rawget(L, 1);
  PCcontinuator* pcc = (PCcontinuator*)lua_topointer(L, -1);
  lua_pop(L, 2);
  
  pcc->set_stepsize(_ssize);
  if(lua_gettop(L) != 0) {
    lua_pushstring(L, "lua_continuator_set_stepsize exit");
    lua_error(L);
  }
  return 0;
}

int lua_continuator_get_stepsize(lua_State* L)
{
  if(lua_gettop(L) != 1) {
    lua_pushstring(L, "lua_continuator_get_stepsize entry");
    lua_error(L);
  }
  lua_pushstring(L, "ud");
  lua_rawget(L, 1);
  PCcontinuator* pcc = (PCcontinuator*)lua_topointer(L, -1);
  lua_pop(L, 2);

  lua_pushnumber(L, pcc->get_stepsize());
  
  if(lua_gettop(L) != 1) {
    lua_pushstring(L, "lua_continuator_get_stepsize exit");
    lua_error(L);
  }
  
  return 1;
}

int lua_continuator_progress(lua_State* L)
{
  const int NARG = lua_gettop(L);
  if(NARG < 4) {
    lua_pushstring(L, "lua_continuator_progress entry");
    lua_error(L);
  }

  lua_pushstring(L, "ud");
  lua_rawget(L, 1);
  PCcontinuator* pcc = (PCcontinuator*)lua_topointer(L, -1);
  lua_pop(L, 1);

  lua_pushstring(L, "residual");
  lua_rawget(L, 1);
  
  vector x = veclua_tovector(L, 2);
  vector t = veclua_tovector(L, 3);
  nc_dirn dir = (nc_dirn)lua_tonumber(L, 4);
  
  vector o(x.length());
  
  lua_remove(L, 4);
  lua_remove(L, 3);
  lua_remove(L, 2);
  lua_remove(L, 1);
  
  lua_insert(L, 1);
  
  prog_result pr = pcc->progress(x, o, t, (void*)L, dir);
  
  if(NARG-3 != lua_gettop(L)) {
    lua_pushstring(L, "lua_continuator_progress error");
    lua_error(L);
  }
  lua_pop(L, lua_gettop(L));
  
  veclua_pushtable(L, o);
  veclua_pushtable(L, t);
  lua_pushnumber(L, pr);
  
  if(lua_gettop(L) != 3) {
    lua_pushstring(L, "lua_continuator_progress exit");
    lua_error(L);
  }
  
  return 3;
  
}


int lua_continuator_newton(lua_State* L)
{
  const int NARG = lua_gettop(L);
  if(NARG < 2) {
    lua_pushstring(L, "lua_continuator_newton entry");
    lua_error(L);
  }
  printf("Lua continuation newton entry point\n");

  lua_pushstring(L, "ud");
  lua_rawget(L, 1);
  PCcontinuator* pcc = (PCcontinuator*)lua_topointer(L, -1);
  lua_pop(L, 1);

  lua_pushstring(L, "residual");
  lua_rawget(L, 1);
  
  vector x = veclua_tovector(L, 2);
  
  vector o(x.length());
  
  lua_remove(L, 2);
  lua_remove(L, 1);
  
  lua_insert(L, 1);
  
  printf("Lua continuation pcc newton call\n");
  prog_result pr = pcc->newton(x, o, (void*)L);
  
  if(NARG-1 != lua_gettop(L)) {
    lua_pushstring(L, "lua_continuator_newton error");
    lua_error(L);
  }
  lua_pop(L, lua_gettop(L));
  
  veclua_pushtable(L, o);
  lua_pushnumber(L, pr);
  
  if(lua_gettop(L) != 2) {
    lua_pushstring(L, "lua_continuator_newton exit");
    lua_error(L);
  }
  
  return 2;
  
}

int lua_continuator__newindex(lua_State* L)
{
  static const int NARG = 3;
  if(lua_gettop(L) != NARG) {
    lua_pushstring(L, "lua_continuator__newindex entry");
    lua_error(L);
  }
  lua_pop(L, NARG);
  
  return 0;
}


int lua_continuator__index(lua_State* L)
{
  static const int NARG = 2;
  if(lua_gettop(L) != NARG) {
    lua_pushstring(L, "lua_continuator__index entry");
    lua_error(L);
  }
  lua_pop(L, NARG);
  
  return 0;
}

int lua_continuator__gc(lua_State* L)
{
  static const int NARG = 1;
  if(lua_gettop(L) != NARG) {
    lua_pushstring(L, "lua_continuator__gc entry");
    lua_error(L);
  }
  PCcontinuator* garbage = (PCcontinuator*)lua_topointer(L, -1);
  
  garbage->~PCcontinuator();
  
  lua_pop(L, NARG);
  
  return 0;
}
