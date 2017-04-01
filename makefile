# Makefile for freesurface continuator project
#
# Builds the lua libraries required for free surface continuation
# with integral methods based on GSL special functions
#
# Revision notes:
#
# 1/3/2013 Phased out wnl and wnl_global, phasing out integrals_hfall, and integrals_triangle.
#
# Phasing out wnl and wnl_global, I believe this stuff will be done with MATLAB
# from now on.

linux : LUA_NAME := lua5.2
macosx : LUA_NAME := lua

LUA_VERSION_REQD = 5.2

LUA_VERSION = $(shell echo "print(_VERSION)" | lua | grep -m 1 -o "[0-9]\.[0-9]" | awk '{print $$1}')

LUA_TEST = $(shell echo $(LUA_VERSION) $(LUA_VERSION_REQD) | awk '{if (($$1 - $$2) >= 0) { print 0 } else { print 1 }}')

linux : LT_OEXT := .so.0.0.0
linux : LT_EXT := .la
linux : LT_PRE := lib
linux : LT_CEXT := .lo
linux : LT_OPRE := .libs/lib
macosx : LT_OEXT := .so
macosx : LT_EXT := .so
macosx : LT_PRE :=
macosx : LT_OPRE :=
macosx : LT_ODIR := ./
macosx : LT_CEXT := .o


#LUA_MODS := integrals_hfall.so integrals_triangle.so \
#  integrals_sbox_homotopy.so \
#  integrals_box_homotopy.so  integrals_box_cluster.so \
#  integrals_pressure_homotopy.so \
#  continuator.so

LUA_MODS := integrals_sbox_homotopy.so \
  integrals_pressure_homotopy.so \
  continuator.so

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
CFLAGS = $(LUA_CFLAGS) $(GSL_CFLAGS) -D __VECTOR_USE_GSL

.PHONY: clean linux macosx checks

default :
	@echo Choose make linux, macosx.


linux : checks memtest $(LUA_MODS)
#	doxygen 2>html/doxygen.wlog > html/doxygen.log

macosx : checks $(LUA_MODS)
#	doxygen 2>html/doxygen.wlog > html/doxygen.log

memtest : memtest.c
	$(CXX) $(CFLAGS) -c -o memtest.o memtest.c
	$(CXX) -o memtest memtest.o $(LUA_LIBS)

checks :
	@echo $(LUA_VERSION)
	@echo $(subst XX, $(LUA_VERSION_REQD), \
  $(subst 0, Lua version XX or greater found, \
  $(subst 1, warning : Require lua version XX or greater, $(LUA_TEST))))

#integrals_hfall.so : integrals_hfall.cpp vector.o
#	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_hfall.cpp
#	$(LUA_MOD_LINK) $(GSL_LIBS) -o $(LT_PRE)integrals_hfall$(LT_EXT) integrals_hfall$(LT_CEXT) vector$(LT_CEXT) 
#	mv $(LT_OPRE)integrals_hfall$(LT_OEXT) integrals_hfall.so

#integrals_triangle.so : integrals_triangle.cpp vector.o
#	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_triangle.cpp
#	$(LUA_MOD_LINK) $(GSL_LIBS) -o $(LT_PRE)integrals_triangle$(LT_EXT) integrals_triangle$(LT_CEXT) vector$(LT_CEXT) 
#	mv $(LT_OPRE)integrals_triangle$(LT_OEXT) integrals_triangle.so

integrals_pressure_homotopy.so : integrals_pressure.o integrals_pressure_formulation.o vector.o integrals_pressure_homotopy.cpp 
	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_pressure_homotopy.cpp
	$(LUA_MOD_LINK) -o $(LT_PRE)integrals_pressure_homotopy$(LT_EXT) integrals_pressure$(LT_CEXT) integrals_pressure_formulation$(LT_CEXT) integrals_pressure_homotopy$(LT_CEXT) vector$(LT_CEXT) $(GSL_LIBS)
	mv $(LT_OPRE)integrals_pressure_homotopy$(LT_OEXT) integrals_pressure_homotopy.so

integrals_pressure_formulation.o : integrals_pressure_formulation.cpp integrals_pressure_formulation.h vector.h
	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_pressure_formulation.cpp

integrals_pressure.o : integrals_pressure.cpp integrals_pressure.h vector.h
	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_pressure.cpp

#integrals_box_cluster.so : integrals_box.o integrals_box_formulation.o vector.o integrals_box_cluster.cpp 
#	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_box_cluster.cpp
#	$(LUA_MOD_LINK) -o $(LT_PRE)integrals_box_cluster$(LT_EXT) integrals_box$(LT_CEXT) integrals_box_formulation$(LT_CEXT) integrals_box_cluster$(LT_CEXT) vector$(LT_CEXT) $(GSL_LIBS)
#	mv $(LT_OPRE)integrals_box_cluster$(LT_OEXT) integrals_box_cluster.so

#integrals_box_homotopy.so : integrals_box.o integrals_box_formulation.o vector.o integrals_box_homotopy.cpp 
#	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_box_homotopy.cpp
#	$(LUA_MOD_LINK) -o $(LT_PRE)integrals_box_homotopy$(LT_EXT) integrals_box$(LT_CEXT) integrals_box_formulation$(LT_CEXT) integrals_box_homotopy$(LT_CEXT) vector$(LT_CEXT) $(GSL_LIBS)
#	mv $(LT_OPRE)integrals_box_homotopy$(LT_OEXT) integrals_box_homotopy.so

#integrals_box_formulation.o : integrals_box_formulation.cpp integrals_box_formulation.h vector.h
#	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_box_formulation.cpp

#integrals_box.o : integrals_box.cpp integrals_box.h vector.h
#	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_box.cpp

integrals_sbox_homotopy.so : integrals_sbox.o integrals_sbox_formulation.o vector.o integrals_sbox_homotopy.cpp 
	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_sbox_homotopy.cpp
	$(LUA_MOD_LINK) -o $(LT_PRE)integrals_sbox_homotopy$(LT_EXT) integrals_sbox$(LT_CEXT) $(GSL_LIBS) integrals_sbox_formulation$(LT_CEXT) integrals_sbox_homotopy$(LT_CEXT) vector$(LT_CEXT) 
	mv $(LT_OPRE)integrals_sbox_homotopy$(LT_OEXT) integrals_sbox_homotopy.so

integrals_sbox_formulation.o : integrals_sbox_formulation.cpp integrals_sbox_formulation.h vector.h
	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_sbox_formulation.cpp

integrals_sbox.o : integrals_sbox.cpp integrals_sbox.h vector.h
	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_sbox.cpp

vector.o : vector.cpp
	$(LUA_MOD_COMPILE) $(CFLAGS) vector.cpp

continuator.so : continuator.cpp vector.o
	$(LUA_MOD_COMPILE) $(CFLAGS) continuator.cpp
	$(LUA_MOD_LINK) -o $(LT_PRE)continuator$(LT_EXT) continuator$(LT_CEXT) vector$(LT_CEXT) $(GSL_LIBS)
	mv $(LT_OPRE)continuator$(LT_OEXT) continuator.so

clean :
	rm -f html/*
	rm -f .libs/*
	rm -f *~
	rm -f *.so *.lo *.la *.o
