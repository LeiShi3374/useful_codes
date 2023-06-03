#!
#! Copyright (c) Stanford University, The Regents of the University of
#!               California, and others.
#!
#! All Rights Reserved.
#!
#! See Copyright-SimVascular.txt for additional details.
#!
#! Permission is hereby granted, free of charge, to any person obtaining
#! a copy of this software and associated documentation files (the
#! "Software"), to deal in the Software without restriction, including
#! without limitation the rights to use, copy, modify, merge, publish,
#! distribute, sublicense, and/or sell copies of the Software, and to
#! permit persons to whom the Software is furnished to do so, subject
#! to the following conditions:
#!
#! The above copyright notice and this permission notice shall be included
#! in all copies or substantial portions of the Software.
#!
#! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
#! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
#! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
#! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
#! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#!
#!--------------------------------------------------------------------

.SUFFIXES: .f .o

CUR_DIR:= $(dir $(abspath $(lastword $(MAKEFILE_LIST)) ))
SRC_DIR = $(CUR_DIR)src
OBJ_DIR = $(CUR_DIR)obj
BIN_DIR = $(CUR_DIR)bin

FC = gfortran
FFLAGS  = -O3 -cpp
FFLAGS += -J $(OBJ_DIR)
LAPACK_LFLAGS = -llapack

debug = 0
ifeq ($(debug),1)
   FFLAGS += -g -fbacktrace -fbounds-check -Wall
endif

S1S2 = 0
ifeq ($(S1S2),1)
   FFLAGS += -DS1S2REST
endif

INCLUDES += -I./include

MYFUN  = COM_TYPE.f    \
         COM_UTIL.f    \
         COM_CHNL.f    \
         COM_LIST.f    \
         MATFUN.f      \
         CEPMOD_AP.f   \
         CEPMOD_BO.f   \
         CEPMOD_FN.f   \
         CEPMOD_TTP.f  \
         ECMOD_DCPLD.f \
         CEPMOD.f      \
         COM_MOD.f     \
         CEPION.f      \
         MAIN.f        \
         READFILES.f

MYEXE  = $(BIN_DIR)/cep.exe

MYSRC = $(patsubst %,   $(SRC_DIR)/%,   $(MYFUN))
MYOBJ = $(patsubst %.f, $(OBJ_DIR)/%.o, $(MYFUN))

OPTIONS = $(FFLAGS)  $(LAPACK_LFLAGS)

#---Rules---#
.PHONY: $(MYEXE)
$(MYEXE): $(MYOBJ)  $(BIN_DIR)
	$(FC) $(MYOBJ) $(OPTIONS) -o $@

$(MYOBJ): | $(OBJ_DIR)

$(OBJ_DIR):
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR)

cleanall:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

