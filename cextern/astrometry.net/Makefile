
all: test_healpix

OBJS = healpix-utils.o healpix.o starutil.o permutedsort.o mathutil.o bl.o qsort_reentrant.o

HEADERS = healpix-utils.h healpix.h

$(OBJS): %.o: %.c $(HEADERS)
	$(CC) -o $@ -c $<

%.o: %.c
	$(CC) -o $@ -c $<

test_healpix: test_healpix-main.c test_healpix.c $(OBJS) CuTest.o

example: example.c $(OBJS)
