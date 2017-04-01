#include "solver.h"

int lua_solver_new(lua_State*);
int lua_solver__gc(lua_State*);
int lua_solver__index(lua_State*);
int lua_solver__newindex(lua_State*);
int lua_solver_progress(lua_State*);
void lua_solver_solvelayer(const vector&, vector&, void*);

struct lua_constant {
  const char *str;
  int       value;
};

static struct lua_constant solver_constants[] = {
  {"CONVERGED", SOLVER_CONVERGED},
  {"DIVERGED", SOLVER_DIVERGED},
  {"ERROR", SOLVER_ERROR},
  {NULL, 0}
};

static const luaL_reg solverR[] = {
  {"new", lua_solver_new},
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

int luaopen_solver(lua_State *L) {
  int NARG = lua_gettop(L);
  luaL_register(L, "solver", solverR);
  
  add_constants(L, solver_constants);
  
  if(lua_gettop(L) != NARG+1) {
    lua_pushstring(L, "luaopen_solver exit");
    lua_error(L);
  }
  return 1;
}

/** Finds near-exact root of an equation \f$H(\mathbf{x}) = \mathbf{0}\f$ by
    Newton iteration.
    
    \param x an (n x 1) vector of the initial approximate solution
    \param o an (n x 1) vector to store the final approximate solution
    \param args extra arguments to pass to H
    \param tol specified accuracy
  
    This function uses Newton's method by solving the following system to find
    the next iteration
    \f[ J(\mathbf{x}_n) \cdot (\mathbf{x}_{n+1} - \mathbf{x}_n) =
        -H(\mathbf{x}_n) \f]
    where \f$J\f$ is the Jacobian matrix of \f$H\f$.
    Jacobian \f$J\f$ used for Newton's method is calculated via the 
    finite difference approximation 
    \f[ J = [J_{ij}] = \frac{\partial H_i}{\partial x_j} \approx
       \frac{H_i(\mathbf{x} + h \mathbf{e}_j) - H_i(\mathbf{x})}{h} \f] where
    \f$H_i(\mathbf{x}) = H(\mathbf{x}) \cdot \mathbf{e}_i\f$

    This method is iterated until either a maximum number of attempts is
    reached (SOLVER_MAX_ATTEMPTS), or the specified accuracy is obtained, that
    is
    \f$ \mathcal{L}_{\mathrm{sup}} H(\mathbf{x}) < tol \f$
*/
solve_result Newtonsolver::converge(const vector& x,
                                    vector&       o,
                                    void*      args,
                                    const real  tol) const
{
  unsigned int      n_col = x.length();
  unsigned int   attempts = 0;
  real            err_val = tol;
  solve_result      r_val = SOLVER_CONVERGED;

  real *heap_v = new real[(n_row*2)+n_col];
  real *heap_m = new real[n_row*n_col];

  vector f_o  = vector(n_row, &(heap_v[0]));
  vector f_x1 = vector(n_row, &(heap_v[n_row]));
  vector x1   = vector(n_col, &(heap_v[n_row*2]));

  matrix J = matrix(n_row, n_col, &(heap_m[0]));

  for(unsigned int i = 1; i <= n_col; i++) o(i) = x(i);

  H(o, f_o, args);

  while(err_val >= tol) {
    if(attempts >= SOLVER_MAX_ATTEMPTS) {
      r_val = SOLVER_NOTCONVERGED;
      break;
    }

    for(unsigned int i = 1; i <= n_col; i++) x1(i) = o(i);

    for(unsigned int j = 1; j <= n_col; j++) {

      x1(j) = o(j) + h;
      H(x1, f_x1, args);
      x1(j) = o(j);

      for(unsigned int i = 1; i <= n_row; i++) J(i,j) = (f_x1(i) - f_o(i)) / h;
    }

    // Solve J dx = -F
    for(unsigned int i = 1; i <= n_row; i++) f_o(i) = -f_o(i);

    solve_lin_mat(J, x1, f_o);

    for(unsigned int i = 1; i <= n_col; i++) o(i) = o(i) + x1(i);
 
    // check if x has become singular
    if(o.issingular()) {
      r_val = SOLVER_DIVERGED;
      break;
    }

    H(o, f_o, args);

    if(f_o.issingular()) {
      r_val = SOLVER_DIVERGED;
      break;
    }

    err_val = 0.0;
    for(unsigned int i = 1; i <= n_col; i++)
      err_val = (err_val < fabs(f_o(i))) ? fabs(f_o(i)) : err_val;
    attempts++;
  }

  delete [] heap_v;
  delete [] heap_m;

  return r_val;

}

/** Finds near-exact root of an equation by Newton iteration.

    \param x the initial approximation
    \param o the final approximate
    \param args extra arguments to pass to H

    This function uses Newton's method by solving the following system to find
    the next iteration
    \f[ J(\mathbf{x}_n) \cdot (\mathbf{x}_{n+1} - \mathbf{x}_n) =
        -H(\mathbf{x}_n) \f]
    where \f$J\f$ is the Jacobian matrix of \f$H\f$.
    Jacobian \f$J\f$ used for Newton's method is calculated via the 
    finite difference approximation 
    \f[ J = [J_{ij}] = \frac{\partial H_i}{\partial x_j} \approx
       \frac{H_i(\mathbf{x} + h \mathbf{e}_j) - H_i(\mathbf{x})}{h} \f] where
    \f$H_i(\mathbf{x}) = H(\mathbf{x}) \cdot \mathbf{e}_i\f$

    This method is iterated until either a maximum number of attempts is
    reached (SOLVER_MAX_ATTEMPTS), or the specified accuracy is obtained, that
    is
    \f$ \mathcal{L}_{\mathrm{sup}}H(\mathbf{x}) < \mathrm{SOLVER_DEFAULT_TOL}\f$
*/
solve_result Newtonsolver::converge(const vector& x,
                                    vector&       o,
                                    void*      args) const
{
  return converge(x, o, args, SOLVER_DEFAULT_TOL);
}

/** Finds near-exact root of an equation by Newton iteration.

    \param x the vector for starting and final approximate
    \param args extra arguments to pass to H

    This function uses Newton's method by solving the following system to find
    the next iteration
    \f[ J(\mathbf{x}_n) \cdot (\mathbf{x}_{n+1} - \mathbf{x}_n) =
        -H(\mathbf{x}_n) \f]
    where \f$J\f$ is the Jacobian matrix of \f$H\f$.
    Jacobian \f$J\f$ used for Newton's method is calculated via the 
    finite difference approximation 
    \f[ J = [J_{ij}] = \frac{\partial H_i}{\partial x_j} \approx
       \frac{H_i(\mathbf{x} + h \mathbf{e}_j) - H_i(\mathbf{x})}{h} \f] where
    \f$H_i(\mathbf{x}) = H(\mathbf{x}) \cdot \mathbf{e}_i\f$

    This method is iterated until either a maximum number of attempts is
    reached (SOLVER_MAX_ATTEMPTS), or the specified accuracy is obtained, that
    is
    \f$ \mathcal{L}_{\mathrm{sup}}H(\mathbf{x}) < \mathrm{SOLVER_DEFAULT_TOL}\f$
*/
solve_result Newtonsolver::converge(vector& x, void* args) const
{
  return converge(x, x, args, SOLVER_DEFAULT_TOL);
}

/** Improves the solution of an approximate root
    \param x the initial approximation
    \param o the improved approximate
    \param args extra arguments to pass to H

    This function finds the next iteration in Newton's method
    \f[ J(\mathbf{x}_n) \cdot (\mathbf{x}_{n+1} - \mathbf{x}_n) =
        -H(\mathbf{x}_n) \f]
    where \f$J\f$ is the Jacobian matrix of \f$H\f$.
    Jacobian \f$J\f$ is calculated via the finite difference approximation 
    \f[ J = [J_{ij}] = \frac{\partial H_i}{\partial x_j} \approx
       \frac{H_i(\mathbf{x} + h \mathbf{e}_j) - H_i(\mathbf{x})}{h} \f] where
    \f$H_i(\mathbf{x}) = H(\mathbf{x}) \cdot \mathbf{e}_i\f$
*/
int Newtonsolver::improve(const vector& x, vector& o, void* args) const
{
  unsigned int n_col = x.length();

  real *heap_v = new real[(n_row*2)+n_col];
  real *heap_m = new real[n_row*n_col];

  vector f_o  = vector(n_row, &(heap_v[0]));
  vector f_x1 = vector(n_row, &(heap_v[n_row]));
  vector x1   = vector(n_col, &(heap_v[n_row*2]));

  matrix J = matrix(n_row, n_col, &(heap_m[0]));

  for(unsigned int i = 1; i <= n_col; i++) {
    o(i) = x(i);
    x1(i) = o(i);
  }

  // Calculate Jacobian matrix

  H(o, f_o, args);

  for(unsigned int j = 1; j <= n_col; j++) {

    x1(j) = o(j) + h;
    H(x1, f_x1, args);
    x1(j) = o(j);

    for(unsigned int i = 1; i <= n_row; i++) J(i,j) = (f_x1(i) - f_o(i)) / h;

  }

  // Solve J dx = -F
  for(unsigned int i = 1; i <= n_row; i++) f_o(i) = -f_o(i);

  solve_lin_mat(J, x1, f_o);

  for(unsigned int i = 1; i <= n_col; i++) o(i) = o(i) + x1(i);

  delete [] heap_v;
  delete [] heap_m;

  return 0;
}

int lua_solver_new(lua_State* L)
{
  static const int NARG = 2;
  int n_rows;
  
  if(lua_gettop(L) != NARG) {
    lua_pushstring(L, "lua_solver_new entry");
    lua_error(L);
  }
  
  switch(lua_type(L, 1))
  {
  case LUA_TFUNCTION :
    lua_newtable(L);

    lua_pushstring(L, "residual");
    lua_pushvalue(L, 1);
    lua_rawset(L, -3);
    
    lua_remove(L, 1);
        
    break;
  default:
    lua_pop(L, NARG);
    return 0;
  }
  
  switch(lua_type(L, 1))
  {
  case LUA_TNUMBER :
    n_rows = lua_tonumber(L, 1);
    lua_remove(L, 1);
        
    break;
  default:
    lua_pop(L, NARG);
    return 0;
  }
  
  lua_pushstring(L, "solve");
  lua_pushcfunction(L, lua_solver_solve);
  lua_rawset(L, -3);

  lua_pushstring(L, "ud");

  // make userdata w/- custom metatable

  void* buffer = lua_newuserdata(L, sizeof(Newtonsolver));
  
  Newtonsolver* nc = new(buffer) Newtonsolver(lua_solver_solvelayer, n_rows);
  
  /** \todo metatable for the userdata */
  lua_newtable(L);
  
  lua_pushstring(L, "__gc");
  lua_pushcfunction(L, lua_solver__gc);
  lua_rawset(L, -3);

  lua_pushstring(L, "__newindex");
  lua_pushcfunction(L, lua_solver__newindex);
  lua_rawset(L, -3);

  lua_pushstring(L, "__index");
  lua_pushcfunction(L, lua_solver__index);
  lua_rawset(L, -3);

  lua_pushstring(L, "__metatable");
  lua_pushvalue(L, -2);
  lua_rawset(L, -3);
  
  lua_setmetatable(L, -2);
  
  lua_rawset(L, -3);
  
  if(lua_gettop(L) != 1) {
    lua_pushstring(L, "lua_solver_new exit");
    lua_error(L);
  }
  
  return 1;
}

void lua_solve_solvelayer(const vector& input, vector& output, void* args)
{
  lua_State* L = (lua_State*)args;
  
  const int NARG = lua_gettop(L);
  
  lua_pushvalue(L, 1);
  veclua_pushtable(L, input);
  for(unsigned int i = 2; i <= NARG; i++) {
    lua_pushvalue(L, i);
  }
  lua_call(L, NARG, 1);

  output = veclua_tovector(L, -1);
  lua_pop(L, 1);
  
  if(NARG != lua_gettop(L)) {
    lua_pushstring(L, "lua_solve_solvelayer exit");
    lua_error(L);
  }
}


int lua_solver_solve(lua_State* L)
{
  const int NARG = lua_gettop(L);
  if(NARG < 3) {
    lua_pushstring(L, "lua_solver_solve exit");
    lua_error(L);
  }

  lua_pushstring(L, "ud");
  lua_rawget(L, 1);
  Newtonsolver* ns = (ns*)lua_topointer(L, -1);
  lua_pop(L, 1);

  lua_pushstring(L, "residual");
  lua_rawget(L, 1);
  
  vector x = veclua_tovector(L, 2);
  nc_dirn dir = (nc_dirn)lua_tonumber(L, 3);
  
  vector o(x.length());
  
  lua_remove(L, 3);
  lua_remove(L, 2);
  lua_remove(L, 1);
  
  lua_insert(L, 1);
  
  solver_result sr = sr->converge(x, o, (void*)L, dir);
  
  if(NARG-2 != lua_gettop(L)) {
    lua_pushstring(L, "lua_solver_solve error");
    lua_error(L);
  }
  lua_pop(L, lua_gettop(L));
  
  veclua_pushtable(L, o);
  lua_pushnumber(L, pr);
  
  if(lua_gettop(L) != 2) {
    lua_pushstring(L, "lua_solver_solve exit");
    lua_error(L);
  }
  
  return 2;
  
}

int lua_solver__newindex(lua_State* L)
{
  static const int NARG = 3;
  if(lua_gettop(L) != NARG) {
    lua_pushstring(L, "lua_solver__newindex entry");
    lua_error(L);
  }
  lua_pop(L, NARG);
  
  return 0;
}

int lua_solver__index(lua_State* L)
{
  static const int NARG = 2;
  if(lua_gettop(L) != NARG) {
    lua_pushstring(L, "lua_solver__index entry");
    lua_error(L);
  }
  lua_pop(L, NARG);
  
  return 0;
}

int lua_solver__gc(lua_State* L)
{
  static const int NARG = 1;
  if(lua_gettop(L) != NARG) {
    lua_pushstring(L, "lua_solver__gc entry");
    lua_error(L);
  }
  Newtonsolver* garbage = (Newtonsolver*)lua_topointer(L, -1);
  
  garbage->~Newtonsolver();
  
  lua_pop(L, NARG);
  
  return 0;
}
