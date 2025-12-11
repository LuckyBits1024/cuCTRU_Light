CC ?= /usr/bin/cc
CFLAGS += -Wall -march=native -mtune=native -O3 -fomit-frame-pointer -Wno-unknown-pragmas
LDFLAGS += -lcrypto -ldl

SOURCES = ctru.c poly.c reduce.c coding.c randombytes.c cbd.c ntt.c inverse.c\
  cpucycles.c speed.c pack.c kem.c fips202.c
HEADERS = params.h ctru.h poly.h reduce.h coding.h ntt.h inverse.h\
  randombytes.h cpucycles.h speed.h pack.h api.h cbd.h fips202.h crypto_hash_sha3256.h symmetric_crypto.h
SOURCES_PRINT = print_ctru.c print_kem.c poly.c reduce.c coding.c randombytes.c cbd.c ntt.c inverse.c\
  cpucycles.c speed.c pack.c fips202.c

all: test_ntt test_pkelight test_kemlight test_speedlight test_printlight


test_ntt: test_ntt.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS)  $< $(SOURCES) -o $@ $(LDFLAGS)

test_pkelight: test_pke.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS)  $< $(SOURCES) -o $@ $(LDFLAGS)

test_kemlight: test_kem.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) $< $(SOURCES) -o $@ $(LDFLAGS)

test_speedlight: test_speed.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) $< $(SOURCES) -o $@ $(LDFLAGS)

test_printlight: test_print.c $(SOURCES_PRINT) $(HEADERS) 
	$(CC) $(CFLAGS) $< $(SOURCES_PRINT) -o $@ $(LDFLAGS)

.PHONY: clean

clean:
	rm -f *~ *~
	rm -f test_kemlight
	rm -f test_speedlight	
	rm -f test_ntt
	rm -f test_pkelight
	rm -f test_printlight

