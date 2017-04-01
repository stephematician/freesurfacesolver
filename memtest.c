#include "stdio.h"

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}

int main(int argc, char **argv)
{
 
  lua_State *L = lua_open();
  printf("lua opened succesfully\n");

  luaL_openlibs(L);
  printf("lua libs opened\n");

  if(argc > 1) {
    // Expect the first argument to be the lua script to run

    int rval = luaL_loadfile(L, argv[1]);
    
    if(!rval) {
      rval |= lua_pcall(L, 0, LUA_MULTRET, 0);
    }
    
    if(!rval) {
      printf("lua script ran successfully\n");
    } else {
      printf("lua script failed\n");
    }
    
  }
  
  lua_close(L);

  printf("lua closed successfully\n");

  return 0; 
}
