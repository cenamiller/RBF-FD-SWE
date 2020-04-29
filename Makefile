

TARGET = pdex

ifeq ($(COMPILER),gnu)
	CC = gcc 
	CFLAGS = -std=c99 -O3 -ffast-math -Wall
	INCLUDES = -I include
	LINKFLAGS = -Wall
	LIBS = -lm
else ifeq ($(COMPILER),qcc)
        CC = /opt/cacheq/qcc/0.5/bin/qcc
        CFLAGS = -e RK4_advance 
        INCLUDES = -I include
        LINKFLAGS = -Wall
        LIBS = -lm
else
	CC = icc
	CFLAGS = -std=c99 -O3 -fno-alias -xHost 
	INCLUDES = -I.
	LINKFLAGS =
	LIBS = -lm 
endif

.PHONY: clean all default

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) $(LINKFLAGS) $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f *.cq
	-rf -f *.gcc
	-rm -f $(TARGET)
