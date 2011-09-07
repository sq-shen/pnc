CXX := g++
SYMBOLS :=
CFLAGS :=
INCLUDES := -I../../comm
LDFLAGS :=
LIBS := -lm -litpp -lfftw3 -llapack -lblas

DEBUG = 1

ifeq ($(DEBUG),1)
	CFLAGS  += -g -Wall -DNDEBUG -pipe
	
	SYMBOLS += -DDBGSHOW
	SYMBOLS += -DCHECKNUM
	SYMBOLS += -DCOMM_ASSERT
	SYMBOLS += -DDBGFILE
endif

ifeq ($(DEBUG),0)
	CFLAGS  += -Wall -O3
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

PNCML_PROG := pnc_ml
PNCML_OBJS := pnc_ml.o

PNCMMSE_PROG := pnc_mmse
PNCMMSE_OBJS := pnc_mmse.o

MIMOPNC_PROG := mimo_pnc
MIMOPNC_OBJS := mimo_pnc.o

DRAWCONST_PROG := draw_const
DRAWCONST_OBJS := draw_const.o

GENCH_PROG := gen_ch
GENCH_OBJS := gen_ch.o 

PNCHARD_PROG := pnc_hard
PNCHARD_OBJS := pnc_hard.o

MRCPNC_PROG := mrc_pnc
MRCPNC_OBJS := mrc_pnc.o

TESTEIGEN_PROG := test_eigen
TESTEIGEN_OBJS := test_eigen.o mimomac.o zftwrelay.o

COMM_OBJS += mimomac.o
COMM_OBJS += zftwrelay.o


#all: $(PROGRAM) $(PNCML_PROG) $(PNCMMSE_PROG) $(MIMOPNC_PROG) $(GENCH_PROG) $(PNCHARD_PROG) $(TESTEIGEN_PROG) $(DRAWCONST_PROG) $(MRCPNC_PROG)
all: $(MRCPNC_PROG)

$(BIN_DIR):
	@[ -d $@ ] || mkdir $@ 

#
# pnc_zf
#
$(PROGRAM): $(BIN_DIR) $(BIN_DIR)/$(PROGRAM)

$(BIN_DIR)/$(PROGRAM): $(addprefix $(BIN_DIR)/, $(OBJS)) $(addprefix $(BIN_DIR)/, $(COMM_OBJS))
	$(CXX) $(LIBS) -o $@ $^
	
#
# pnc_ml
#
$(PNCML_PROG): $(BIN_DIR) $(BIN_DIR)/$(PNCML_PROG)

$(BIN_DIR)/$(PNCML_PROG): $(addprefix $(BIN_DIR)/, $(PNCML_OBJS)) $(addprefix $(BIN_DIR)/, $(COMM_OBJS))
	$(CXX) $(LIBS) -o $@ $^
	
#
# pnc_mmse
#
$(PNCMMSE_PROG): $(BIN_DIR) $(BIN_DIR)/$(PNCMMSE_PROG)

$(BIN_DIR)/$(PNCMMSE_PROG): $(addprefix $(BIN_DIR)/, $(PNCMMSE_OBJS)) $(addprefix $(BIN_DIR)/, $(COMM_OBJS))
	$(CXX) $(LIBS) -o $@ $^
	
	
#
# mimo_pnc
#
$(MIMOPNC_PROG): $(BIN_DIR) $(BIN_DIR)/$(MIMOPNC_PROG)

$(BIN_DIR)/$(MIMOPNC_PROG): $(addprefix $(BIN_DIR)/, $(MIMOPNC_OBJS)) $(addprefix $(BIN_DIR)/, $(COMM_OBJS))
	$(CXX) $(LIBS) -o $@ $^
	

#
# mrc_pnc
#
$(MRCPNC_PROG): $(BIN_DIR) $(BIN_DIR)/$(MRCPNC_PROG)

$(BIN_DIR)/$(MRCPNC_PROG): $(addprefix $(BIN_DIR)/, $(MRCPNC_OBJS)) $(addprefix $(BIN_DIR)/, $(COMM_OBJS))
	$(CXX) $(LIBS) -o $@ $^
	
#
# draw_const
#
$(DRAWCONST_PROG): $(BIN_DIR) $(BIN_DIR)/$(DRAWCONST_PROG)

$(BIN_DIR)/$(DRAWCONST_PROG): $(addprefix $(BIN_DIR)/, $(DRAWCONST_OBJS)) $(addprefix $(BIN_DIR)/, $(COMM_OBJS))
	$(CXX) $(LIBS) -o $@ $^


#
# gen_ch
#
$(GENCH_PROG): $(BIN_DIR) $(BIN_DIR)/$(GENCH_PROG)

$(BIN_DIR)/$(GENCH_PROG): $(addprefix $(BIN_DIR)/, $(GENCH_OBJS))
	$(CXX) $(LIBS) -o $@ $^
	

#
# pnc_hard
#
$(PNCHARD_PROG): $(BIN_DIR) $(BIN_DIR)/$(PNCHARD_PROG)

$(BIN_DIR)/$(PNCHARD_PROG): $(addprefix $(BIN_DIR)/, $(PNCHARD_OBJS)) $(addprefix $(BIN_DIR)/, $(COMM_OBJS))
	$(CXX) $(LIBS) -o $@ $^
	

#
# test
#
$(TESTEIGEN_PROG): $(BIN_DIR) $(BIN_DIR)/$(TESTEIGEN_PROG)

$(BIN_DIR)/$(TESTEIGEN_PROG): $(addprefix $(BIN_DIR)/, $(TESTEIGEN_OBJS))
	$(CXX) $(LIBS) -o $@ $^


#======================================
$(BIN_DIR)/%.o: %.cpp
	$(CXX) $(CFLAGS) $(SYMBOLS) $(INCLUDES) -c -o $@ $<


############################################################
# 
############################################################

.PHONY: all clean $(PROGRAM) $(GENCH_PROG) $(PNCHARD_PROG) $(TESTEIGEN_PROG) $(PNCML_PROG) $(PNCMMSE_PROG) $(DRAWCONST_PROG) $(MIMOPNC_PROG) $(MRCPNC_PROG)

	
#clean:
#	@rm -rf dbglog1.txt dbglog2.txt $(BIN_DIR)
clean:
	@rm -rf dbglog1.txt dbglog2.txt $(BIN_DIR)/*.o mrc_pnc
