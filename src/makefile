LUA_MODS := \
    $(LIB_DIR)/integrals_sbox_homotopy.so \
    $(LIB_DIR)/integrals_pressure_homotopy.so \
    $(LIB_DIR)/continuator.so

linux: $(LUA_MODS)
	@echo $(LUA_MODS)

$(LIB_DIR)/integrals_pressure_homotopy.so : \
    $(OBJ_DIR)/integrals_pressure$(LT_CEXT) \
    $(OBJ_DIR)/integrals_pressure_formulation$(LT_CEXT) \
    $(OBJ_DIR)/vector$(LT_CEXT) \
    integrals_pressure_homotopy.cpp 
	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_pressure_homotopy.cpp \
	    -o $(OBJ_DIR)/integrals_pressure_homotopy$(LT_CEXT)
	$(LUA_MOD_LINK) -o $(LIB_DIR)/$(LT_PRE)integrals_pressure_homotopy$(LT_EXT) \
	    $(OBJ_DIR)/integrals_pressure$(LT_CEXT) \
	    $(OBJ_DIR)/integrals_pressure_formulation$(LT_CEXT) \
	    $(OBJ_DIR)/integrals_pressure_homotopy$(LT_CEXT) \
	    $(OBJ_DIR)/vector$(LT_CEXT) \
	    $(GSL_LIBS)
	mv $(LIB_DIR)/$(LT_OPRE)integrals_pressure_homotopy$(LT_OEXT) $@

$(OBJ_DIR)/integrals_pressure_formulation$(LT_CEXT) : \
    integrals_pressure_formulation.cpp \
    $(INCLUDE_DIR)/integrals_pressure_formulation.h \
    $(INCLUDE_DIR)/vector.h
	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_pressure_formulation.cpp -o $@

$(OBJ_DIR)/integrals_pressure$(LT_CEXT) : \
    integrals_pressure.cpp \
    $(INCLUDE_DIR)/integrals_pressure.h \
    $(INCLUDE_DIR)/vector.h
	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_pressure.cpp -o $@

$(LIB_DIR)/integrals_sbox_homotopy.so : \
    $(OBJ_DIR)/integrals_sbox$(LT_CEXT) \
    $(OBJ_DIR)/integrals_sbox_formulation$(LT_CEXT) \
    $(OBJ_DIR)/vector$(LT_CEXT) \
    integrals_sbox_homotopy.cpp 
	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_sbox_homotopy.cpp \
	    -o $(OBJ_DIR)/integrals_sbox_homotopy$(LT_CEXT)
	$(LUA_MOD_LINK) -o $(LIB_DIR)/$(LT_PRE)integrals_sbox_homotopy$(LT_EXT) \
	    $(OBJ_DIR)/integrals_sbox$(LT_CEXT) \
	    $(GSL_LIBS) \
	    $(OBJ_DIR)/integrals_sbox_formulation$(LT_CEXT) \
	    $(OBJ_DIR)/integrals_sbox_homotopy$(LT_CEXT) \
	    $(OBJ_DIR)/vector$(LT_CEXT)
	mv $(LIB_DIR)/$(LT_OPRE)integrals_sbox_homotopy$(LT_OEXT) $@

$(OBJ_DIR)/integrals_sbox_formulation$(LT_CEXT) : \
    integrals_sbox_formulation.cpp \
    $(INCLUDE_DIR)/integrals_sbox_formulation.h \
    $(INCLUDE_DIR)/vector.h
	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_sbox_formulation.cpp -o $@

$(OBJ_DIR)/integrals_sbox$(LT_CEXT) : \
    integrals_sbox.cpp \
    $(INCLUDE_DIR)/integrals_sbox.h \
    $(INCLUDE_DIR)/vector.h
	$(LUA_MOD_COMPILE) $(CFLAGS) integrals_sbox.cpp -o $@

$(OBJ_DIR)/vector$(LT_CEXT) : vector.cpp
	$(LUA_MOD_COMPILE) $(CFLAGS) vector.cpp -o $@

$(LIB_DIR)/continuator.so : continuator.cpp $(OBJ_DIR)/vector$(LT_CEXT)
	$(LUA_MOD_COMPILE) $(CFLAGS) continuator.cpp \
	    -o $(OBJ_DIR)/continuator$(LT_CEXT)
	$(LUA_MOD_LINK) -o $(LIB_DIR)/$(LT_PRE)continuator$(LT_EXT) \
	    $(OBJ_DIR)/continuator$(LT_CEXT) \
	    $(OBJ_DIR)/vector$(LT_CEXT) \
	    $(GSL_LIBS)
	mv $(LIB_DIR)/$(LT_OPRE)continuator$(LT_OEXT) $@

$(BIN_DIR)/memtest : memtest.c
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/memtest.o memtest.c
	$(CXX) -o $(BIN_DIR)/memtest $(OBJ_DIR)/memtest.o $(LUA_LIBS)


clean:
	rm -rf $(LIB_DIR)
	rm -rf $(BIN_DIR)
	rm -rf $(OBJ_DIR)

