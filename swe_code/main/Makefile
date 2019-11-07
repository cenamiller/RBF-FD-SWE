include ../../include.mk

MAIN_OBJS = profiling.o runtime_params.o rk4_rbffd_swe.o main.o

MAIN_LIB = swe_main.a

all: $(MAIN_LIB)

$(MAIN_LIB): $(MAIN_OBJS)
	$(LD_R) $^ -o $@
	@$(MAKE) clean

profiling.o: profiling.c
	$(CC) $(CFLAGS) $(INCLUDE_MAIN) $(INCLUDE_COMMON) -c $< -o $@	

runtime_params.o: runtime_params.c
	$(CC) $(CFLAGS) $(INCLUDE_MAIN) $(INCLUDE_COMMON) -c $< -o $@

rk4_rbffd_swe.o: rk4_rbffd_swe.c
	$(CC) $(CFLAGS) $(OPTRPT) $(INCLUDE_MAIN) $(INCLUDE_MPI) $(INCLUDE_COMMON) -c $< -o $@	

main.o: main.c
	$(CC) $(CFLAGS) $(INCLUDE_ALL) -c $< -o $@	

.PHONY: clean
clean:
	rm -f *.o
