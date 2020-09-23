CC = /usr/bin/cc
#CFLAGS = -Wall -g -Wextra -march=native -mtune=native -O3 -fomit-frame-pointer -msse4.1 
CFLAGS = -fPIC -fvisibility=hidden -Werror -Wall -Wextra -Wpedantic -Wstrict-prototypes -Wshadow -Wformat=2 -Wfloat-equal -Wwrite-strings -O3 -fomit-frame-pointer -fdata-sections -ffunction-sections -Wl,--gc-sections -std=gnu11 -MD -MT
NISTFLAGS = -O3 -fomit-frame-pointer -fPIC

OBJ = 	aes.o aes256ctr.o Alg.o fips202.o fips202x4.o ntt.o packing.o poly.o polyvec.o randombytes.o rng.o sign.o
EX_OBJ = KeccakP-1600-times4-SIMD256.o			
SOURCES = 	main.c \
			Alg.c \
			fips202.c \
			randombytes.c \
			Polynomial.c \
			io.c \
			Rounding.c \

HEADERS = 	randombytes.h \
			fips202.h \
			Parameters.h \
			Alg.h \
			Polynomial.h \
			io.h \
			Rounding.h \


all: test libaigis.a

libaigis.a: $(OBJ) $(EX_OBJ)
	ar csr $@ $^
	rm -rf $^
$(OBJ): %.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(EX_OBJ): %.o: keccak4x/%.c
	$(CC) $(CFLAGS) -c $< -o $@
	
test: test.c api.h libaigis.a 
	$(CC) $(CFLAGS) test.c -o $@ -laigis -L.

.PHONY: clean

clean:
	rm -f test libaigis.a
