# Makefile created by Joao Francisco B. S. Martins <joaofbsm@dcc.ufmg.br>

CFLAGS = -Wall -g

LIBS = -lm

BIN = simplex

OBJS = lalgebra.o

# This default rule compiles the executable program
$(BIN): $(OBJS) $(BIN).c
	$(CC) $(CFLAGS) $(LIBS) $(OBJS) $(BIN).c -o $(BIN)

%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *~ *.o $(BIN)