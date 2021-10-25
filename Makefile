CC = /opt/build/openmpi-4.0.4/bin/mpicc
LD = /opt/build/openmpi-4.0.4/bin/mpicc

CFLAGS   = -I/opt/build/pari-2.13.2/include -I/opt/build/openmpi-4.0.4/include
LDFLAGS  = -I/opt/build/pari-2.13.2/include -I/opt/build/openmpi-4.0.4/include
LIBS     = -L/opt/build/pari-2.13.2/lib -lpari -L/opt/build/openmpi-4.0.4/lib -lmpi

SDIR = src
ODIR = objs
BDIR = .
TDIR = tests

_BINS = main make-integral-posdef-symm-matrix make-diagonal-matrix print-integral-matrix
BINS = $(patsubst %,$(BDIR)/%,$(_BINS))
_TESTS =
TESTS = $(patsubst %,$(BDIR)/%,$(_TESTS))
_OBJS = common.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

all: $(BINS)

$(BDIR)/%: $(OBJS) $(ODIR)/%.o
	$(LD) -g -Wall -o $@ $^ $(LDFLAGS) $(RUNPTH) $(LIBS)

$(ODIR)/%.o: $(SDIR)/%.c $(SDIR)/%.h
	$(CC) -g -Wall -c $(CFLAGS) $< -o $@

$(ODIR)/%.o: $(TDIR)/%.c $(TDIR)/%.h
	$(CC) -g -Wall -c $(CFLAGS) $< -o $@

clean:
	rm -f $(SDIR)/*~ $(ODIR)/*~ $(BDIR)/*~ $(TDIR)/*~ $(OBJS) $(BINS) $(TESTS)
