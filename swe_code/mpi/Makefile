include ../../include.mk

INCLUDE = -I./include $(INCLUDE_COMMON)

MPI_OBJS = init_patches.o halos.o

MPI_LIB = swe_mpi.a

all: $(MPI_LIB)

$(MPI_LIB): $(MPI_OBJS)
	$(LD_R) $^ -o $@
	@$(MAKE) clean

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o

