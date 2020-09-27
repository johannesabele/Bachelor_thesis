BINARY = HEART
CC = gcc
CFLAGS = -O2 -flto
LFLAGS = -lm 
OBJECTS = main.o initialcond.o vectorcalc.o arrayhandling.o mechanical.o helper.o
HEADER = heart.h

$(BINARY): $(OBJECTS)
	$(CC) $(OBJECTS) $(CFLAGS) $(LFLAGS) -o $(BINARY)

%.o : %.c $(HEADER)
	$(CC) $(CFLAGS) $< -c 

run: HEART
	./HEART