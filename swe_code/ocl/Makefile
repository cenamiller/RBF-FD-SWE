include ../../include.mk

INCLUDE = -I./include $(INCLUDE_COMMON)

OCL_OBJS = device_setup.o RK_ocl.o buffers.o

OCL_LIB = swe_ocl.a

all: $(OCL_LIB)

$(OCL_LIB): $(OCL_OBJS)
	$(LD_R) $^ -o $@
	@$(MAKE) clean

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o

