CC=    gcc
LD=    gcc


ifndef PROJ_DIR
PROJ_DIR = Ia_Classic
endif


CFLAGS= -O3 -mavx2 -std=c11 -Wall -Wextra -fno-omit-frame-pointer
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
