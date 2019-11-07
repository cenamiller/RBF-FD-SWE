include ../../include.mk

INCLUDE	= -I./include $(INCLUDE_COMMON)

RCM_OBJS = rcm.o reorder_nodes.o

RCM_LIB	= swe_rcm.a


all: $(RCM_LIB)

$(RCM_LIB): $(RCM_OBJS)
	$(LD_R) $^ -o $@
	@$(MAKE) clean

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o
