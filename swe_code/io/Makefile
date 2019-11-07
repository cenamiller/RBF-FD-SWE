include ../../include.mk

INCLUDE	= -I./include $(INCLUDE_COMMON) $(INCLUDE_NETCDF)

IO_OBJS = input.o

IO_LIB 	= swe_io.a


all: $(IO_LIB)

$(IO_LIB): $(IO_OBJS)
	$(LD_R) $< -o $@
	@$(MAKE) clean

input.o: input.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@	

nc2bin.exe: nc2bin.c input.o
	$(CC) $(CFLAGS) $(INCLUDE) $(NETCDF_LIBS) -g $^ -o $@


.PHONY: clean cleanobjs
clean:
	rm -f *.o

