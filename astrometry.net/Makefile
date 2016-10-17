
all: test_healpix

OBJS = healpix-utils.o healpix.c

HEADERS = healpix-utils.h healpix.h

$(OBJS): %.o: %.c $(HEADERS)
	$(CC) -o $@ -c $<

%.o: %.c
	$(CC) -o $@ -c $<

test_healpix: test_healpix.c $(OBJS) cutest.o
