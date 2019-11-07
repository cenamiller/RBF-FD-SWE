include ../../include.mk

INCLUDE	= -I./include $(INCLUDE_COMMON)

DL_OBJS = matrix_transformations.o layout.o

DL_LIB	= swe_dl.a


all: $(DL_LIB)

$(DL_LIB): $(DL_OBJS)
	$(LD_R) $^ -o $@
	@$(MAKE) clean

test: $(DL_OBJS) test.c
	$(CC) $(CFLAGS) $(INCLUDE) $^ -o $@.exe

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o *.exe

