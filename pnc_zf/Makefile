CXX := g++
SYMBOLS :=
CFLAGS :=
INCLUDES := -I../../comm
LDFLAGS :=
LIBS := -lm

DEBUG = 1

ifeq ($(DEBUG),1)
	CFLAGS  += -g -Wall
	
	SYMBOLS += -DDBGSHOW
	SYMBOLS += -DCHECKNUM
	SYMBOLS += -DCOMM_ASSERT
	#SYMBOLS += -DDBGFILE
endif

ifeq ($(DEBUG),0)
	CFLAGS  += -Wall
endif



COMM_DIR = ../../comm



VPATH = $(COMM_DIR)

############################################################
# Project: pnc_zf
############################################################
PROGRAM = pnc_zf
BIN_DIR = bin
OBJS = 
COMM_OBJS = 

all: $(BIN_DIR) $(PROGRAM)

$(BIN_DIR):
	@[ -d $@ ] || mkdir $@ 

# Main simulation program
OBJS += main.o

COMM_OBJS += channel.o
COMM_OBJS += debug_log.o
COMM_OBJS += modulator.o
COMM_OBJS += random.o
COMM_OBJS += types.o
COMM_OBJS += CRelay.o



$(PROGRAM): $(BIN_DIR)/$(PROGRAM)

$(BIN_DIR)/$(PROGRAM): $(addprefix $(BIN_DIR)/, $(OBJS)) $(addprefix $(BIN_DIR)/, $(COMM_OBJS))
	$(CXX) $(LIBS) -o $@ $^
	
$(BIN_DIR)/%.o: %.cpp
	$(CXX) $(CFLAGS) $(SYMBOLS) $(INCLUDES) -c -o $@ $<


############################################################
# 
############################################################


.PHONY: all clean $(PROGRAM)
	
clean:
	@rm -rf dbglog1.txt dbglog2.txt $(BIN_DIR)
