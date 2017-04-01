# Makefile for freesurface continuator project
#
# Builds the lua libraries required for free surface continuation
# with integral methods based on GSL special functions
#
# Revision notes:
#
# 1/3/2013 Phased out wnl and wnl_global, phasing out integrals_hfall, and integrals_triangle.
# 1/4/2017 Phased out integrals_hfall, integrals_box and integrals_triangle

linux : LUA_NAME := lua5.2
macosx : LUA_NAME := lua

LUA_VERSION_REQD = 5.2

LUA_VERSION = $(shell echo "print(_VERSION)" | lua | grep -m 1 -o "[0-9]\.[0-9]" | awk '{print $$1}')

LUA_TEST = $(shell echo $(LUA_VERSION) $(LUA_VERSION_REQD) | awk '{if (($$1 - $$2) >= 0) { print 0 } else { print 1 }}')

linux : LT_OEXT := .so.0.0.0
linux : LT_EXT := .la
linux : LT_CEXT := .lo
linux : LT_PRE := lib
linux : LT_OPRE := .libs/lib
macosx : LT_OEXT := .so
macosx : LT_EXT := .so
macosx : LT_CEXT := .o
macosx : LT_PRE :=
macosx : LT_OPRE := 

ROOT_DIR    = $(PWD)
LIB_DIR     = $(ROOT_DIR)/lib
INCLUDE_DIR = $(ROOT_DIR)/include
OBJ_DIR     = $(ROOT_DIR)/build
BIN_DIR     = $(ROOT_DIR)/bin

LUA_CFLAGS = $(shell pkg-config --cflags $(LUA_NAME))
LUA_LIBS = $(shell pkg-config --libs $(LUA_NAME))
GSL_VERSION = $(shell pkg-config --modversion gsl)
GSL_CFLAGS = $(shell pkg-config --cflags gsl)
GSL_LIBS = $(shell pkg-config --libs gsl)
INSTALL_CMOD = $(shell pkg-config --variable=INSTALL_CMOD $(LUA_NAME))

linux : LUA_MOD_LINK = libtool --tag=CXX --mode=link g++ -shared -rpath $(INSTALL_CMOD)
linux : LUA_MOD_COMPILE = libtool --tag=CXX --mode=compile g++ -c -shared
macosx : LUA_MOD_LINK = $(CXX) -bundle -undefined dynamic_lookup
macosx : LUA_MOD_COMPILE = $(CXX) -fPIC -c

LIBS = $(LUA_LIBS) $(GSL_LIBS)
CFLAGS = -I$(INCLUDE_DIR) $(LUA_CFLAGS) $(GSL_CFLAGS) -D __VECTOR_USE_GSL

export LT_OEXT LT_EXT LT_CEXT LT_PRE LT_OPRE
export CFLAGS LIBS INCLUDE_DIR OBJ_DIR LIB_DIR BIN_DIR
export LUA_MOD_COMPILE LUA_MOD_LINK

.PHONY: clean linux macosx checks

default :
	@echo Choose make linux, macosx.

linux : checks
	cd src && $(MAKE) $@
#	doxygen 2>html/doxygen.wlog > html/doxygen.log

macosx : checks
	cd src && $(MAKE) $@
#	doxygen 2>html/doxygen.wlog > html/doxygen.log

checks :
	@echo $(LUA_VERSION)
	@echo $(subst XX, $(LUA_VERSION_REQD), \
  $(subst 0, Lua version XX or greater found, \
  $(subst 1, warning : Require lua version XX or greater, $(LUA_TEST))))

clean :
	cd src/ && $(MAKE) clean
	rm -f html/*
	rm -f .libs/*
	rm -f *.so *.lo *.la *.o
