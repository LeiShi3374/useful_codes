#--------------------------------------------------------------------
#
#	This is the Makefile to build LPN for Norwood circulation
#
#--------------------------------------------------------------------

include Makefile.in

MYFFUN   = MOD_TYPE.f \
           MOD_UTIL.f \
           MOD_CHNL.f \
           MOD_LIST.f \
           MOD_COM.f \
           MOD_MATFUN.f \
           MOD_ALLFUN.f \
           NN.f \
           READFILES.f \
           MAIN.f

MYFCFUN  = vtkComMod.f90          \
           vtkTypeParams.f90      \
           vtkXMLParser.f90

MYCFUN   = vtkZpipe.c

MYEXEC   = $(BIN_DIR)/calc-disps.exe

MYSRC    = $(patsubst %,  $(SRC_DIR)/%,  $(MYFFUN))
MYOBJ    = $(patsubst %.f90,$(OBJ_DIR)/%.o,$(MYFCFUN))
MYOBJ   += $(patsubst %.c,$(OBJ_DIR)/%.o,$(MYCFUN))
MYOBJ   += $(patsubst %.f,$(OBJ_DIR)/%.o,$(MYFFUN))

LFLAGS   = $(FCFLAGS) $(Z_LFLAGS)

INCLUDES = $(INC_DIR)
LIBS     = -llapack
OPTIONS  = $(LIBS) $(FFLAGS) $(LFLAGS)

#############################################################################
#---- Rules ----#

.PHONY: $(MYEXEC)
$(MYEXEC): $(MYOBJ)
	$(FC) $^ $(OPTIONS) -o $@

$(MYOBJ): | $(OBJ_DIR) $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR)

cleanall:
	rm -rf $(BIN_DIR) $(OBJ_DIR)

