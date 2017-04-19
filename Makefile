CC=gcc
CFLAGS=-lm -O2 -funroll-loops
EXEC=Run
$(EXEC): *.c *.h
	$(CC) -o $@ Main.c $(CFLAGS)

debug: *.c *.h
	$(CC) -o $@ Main.c $(CFLAGS) -g

clean:
	rm -f $(EXEC) debug
