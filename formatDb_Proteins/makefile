CC = gcc
CFLAGS = -Wall -g
LFLAGS = -lm
RUTAUTIL = util

all: fdb

fdb: FormatDB.o bitLib.o utils.o
	$(CC) $(CFLAGS) $(LFLAGS) FormatDB.o bitLib.o utils.o -o blast_fdb

util: utils.o

bitLib: bitLib.o

# util
utils.o: $(RUTAUTIL)/utils.c $(RUTAUTIL)/utils.h
	$(CC) $(CFLAGS) -c $(RUTAUTIL)/utils.c

#bitlib
bitLib.o: $(RUTAUTIL)/bitLib.c $(RUTAUTIL)/bitLib.h
	$(CC) $(CFLAGS) -c $(RUTAUTIL)/bitLib.c
	
FormatDB.o: FormatDB.c
	$(CC) $(CFLAGS) -c FormatDB.c

clean:
	rm *.o blast_fdb
