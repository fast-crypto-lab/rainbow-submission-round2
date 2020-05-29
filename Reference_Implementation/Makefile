CC=    gcc
LD=    gcc


#
# The variable `$PROJ_DIR' controls the variant(corresponding to a specific directory) will be built.
# To build a specific variant, set the $PROJ_DIR to a specific name of the directory.
#
# For example of building the `Ia_Classic' variant,
# one have to set $PROJ_DIR = Ia_Classic.
# The makefile will compile codes in the `Ia_Classic' directory only.
# No code changes have to be made.
#
# All possible variants are listed as followings.
#

ifndef PROJ_DIR
PROJ_DIR = Ia_Classic
#PROJ_DIR = Ia_Cyclic
#PROJ_DIR = Ia_CompressedCyclic
#PROJ_DIR = IIIc_Classic
#PROJ_DIR = IIIc_Cyclic
#PROJ_DIR = IIIc_CompressedCyclic
#PROJ_DIR = Vc_Classic
#PROJ_DIR = Vc_Cyclic
#PROJ_DIR = Vc_CompressedCyclic
endif


CFLAGS= -O3 -std=c11 -Wall -Wextra -fno-omit-frame-pointer
INCPATH= -I/usr/local/include -I/opt/local/include -I/usr/include -I$(PROJ_DIR)
LDFLAGS=
LIBPATH= -L/usr/local/lib -L/opt/local/lib -L/usr/lib
LIBS=    -lcrypto


ifeq ($(PROJ_DIR),avx2)
CFLAGS += -mavx2
CXXFLAGS += -mavx2
endif



SRCS = $(wildcard $(PROJ_DIR)/*.c)
SRCS_O = $(SRCS:.c=.o)
SRCS_O_ND = $(subst $(PROJ_DIR)/,,$(SRCS_O))

OBJ = $(SRCS_O_ND)

EXE= rainbow-genkey rainbow-sign rainbow-verify PQCgenKAT_sign


CSRC= $(wildcard *.c)


ifdef DEBUG
        CFLAGS=  -D_DEBUG_ -g -O1 -mavx2 -std=c99 -Wall -Wextra -fsanitize=address -fno-omit-frame-pointer
        CXXFLAGS= -D_DEBUG_ -g -O1 -mavx2 -Wall -Wextra -fno-exceptions -fno-rtti -nostdinc++
endif

ifdef GPROF
	CFLAGS += -pg
	CXXFLAGS += -pg
	LDFLAGS += -pg
endif

.PHONY: all tests tables clean

all: $(OBJ) $(EXE)


%-test: $(OBJ) %-test.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

%-benchmark: $(OBJ) %-benchmark.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

rainbow-genkey: $(OBJ) rainbow-genkey.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

rainbow-sign: $(OBJ) rainbow-sign.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

rainbow-verify: $(OBJ) rainbow-verify.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

PQCgenKAT_sign: $(OBJ) PQCgenKAT_sign.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCPATH) -c $<

%.o: $(PROJ_DIR)/%.c
	$(CC) $(CFLAGS) $(INCPATH) -c $<

clean:
	rm *.o *-test *-benchmark rainbow-genkey rainbow-sign rainbow-verify PQCgenKAT_sign;
