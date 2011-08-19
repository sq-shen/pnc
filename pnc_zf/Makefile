CXX := g++
SYMBOLS :=
CFLAGS :=
INCLUDES := -I../../comm
LDFLAGS :=
LIBS := -lm -litpp -lfftw3 -llapack -lblas

DEBUG = 1

ifeq ($(DEBUG),1)
	CFLAGS  += -g -Wall -DNDEBUG -O3 -pipe
	
	SYMBOLS += -DDBGSHOW
	SYMBOLS += -DCHECKNUM
	SYMBOLS += -DCOMM_ASSERT
	SYMBOLS += -DDBGFILE
endif

ifeq ($(DEBUG),0)
	CFLAGS  += -Wall
endif



COMM_DIR = ../../comm



VPATH = $(COMM_DIR)

############################################################
# Project: pnc_zf
############################################################
BIN_DIR := bin
COMM_OBJS := 

PROGRAM := pnc_zf
OBJS := pnc_zf.o

GENCH_PROG := gen_ch
GENCH_OBJS := gen_ch.o 

COMM_OBJS += debug_log.o
COMM_OBJS += mimomac.o
COMM_OBJS += zftwrelay.o


all: $(PROGRAM) $(GENCH_PROG)

$(BIN_DIR):
	@[ -d $@ ] || mkdir $@ 

#
# pnc_zf
#
$(PROGRAM): $(BIN_DIR) $(BIN_DIR)/$(PROGRAM)

$(BIN_DIR)/$(PROGRAM): $(addprefix $(BIN_DIR)/, $(OBJS)) $(addprefix $(BIN_DIR)/, $(COMM_OBJS))
	$(CXX) $(LIBS) -o $@ $^

#
# gen_ch
#
$(GENCH_PROG): $(BIN_DIR) $(BIN_DIR)/$(GENCH_PROG)

$(BIN_DIR)/$(GENCH_PROG): $(addprefix $(BIN_DIR)/, $(GENCH_OBJS))
	$(CXX) $(LIBS) -o $@ $^


	
$(BIN_DIR)/%.o: %.cpp
	$(CXX) $(CFLAGS) $(SYMBOLS) $(INCLUDES) -c -o $@ $<


############################################################
# 
############################################################

.PHONY: all clean $(PROGRAM) $(GENCH_PROG)
	
clean:
	@rm -rf dbglog1.txt dbglog2.txt $(BIN_DIR)
