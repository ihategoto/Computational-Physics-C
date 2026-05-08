CC      = gcc
CFLAGS  = -Wall -Wextra -O2 -std=c99
INC     = -I include
LIBDIR  = lib
SRCDIR  = src
EXDIR   = examples

SRCS    = $(SRCDIR)/ode.c \
          $(SRCDIR)/ode3d.c \
          $(SRCDIR)/random.c \
          $(SRCDIR)/walker.c \
          $(SRCDIR)/percolation.c
OBJS    = $(SRCS:.c=.o)
LIB     = $(LIBDIR)/libphyslib.a

EXAMPLES = harmonic_oscillator rossler random_walk percolation_ex lattice_gas

.PHONY: all examples clean

all: $(LIB)

$(LIB): $(OBJS) | $(LIBDIR)
	ar rcs $@ $^

$(LIBDIR):
	mkdir -p $(LIBDIR)

$(SRCDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

examples: $(LIB) $(EXAMPLES)

$(EXAMPLES): %: $(EXDIR)/%.c $(LIB)
	$(CC) $(CFLAGS) $(INC) $< -L$(LIBDIR) -lphyslib -lm -o $@

plots: examples
	python3 plots/generate_plots.py

clean:
	rm -f $(SRCDIR)/*.o $(EXAMPLES) $(LIB)
