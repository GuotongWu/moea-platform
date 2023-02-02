# Makefile for Monitor

# Compiler
CC = gcc

# Compiler options
CFLAGS = -g -Wall -pedantic

# all object files
SEL_OBJECTS = monitor_user.o monitor.o monitor_internal.o

monitor : $(SEL_OBJECTS)
	$(CC) $(CFLAGS) -lm $(SEL_OBJECTS) -o monitor

monitor_internal.o : monitor_internal.c monitor_internal.h monitor.h monitor_user.h
	$(CC) $(CFLAGS) -c monitor_internal.c 

monitor_user.o : monitor_user.c monitor_user.h monitor.h
	$(CC) $(CFLAGS) -c monitor_user.c

monitor.o : monitor.c monitor.h monitor_user.h monitor_internal.h
	$(CC) $(CFLAGS) -c monitor.c

clean:
	rm -f *~ *.o
