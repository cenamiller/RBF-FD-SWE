#########################################
#					#
#	Top Makefile for SWE Model	#
#					#
#########################################

include include.mk

all:
	@$(MAKE) clean
	@$(MAKE) swe

debug:
	@$(MAKE) all CFLAGS+=-g

main:
	@$(MAKE) -C $(MAIN_DIR)

rcm:
	@$(MAKE) -C $(RCM_DIR)

layout:
	@$(MAKE) -C $(LAYOUT_DIR)

io:
	@$(MAKE) -C $(IO_DIR)

mpi:
	@$(MAKE) -C $(MPI_DIR)

ocl:
	@$(MAKE) -C $(OCL_DIR)

swe: io rcm layout mpi ocl main
	$(MKDIR_P) $(EXEC_DIR)
	$(CC) $(CFLAGS) $(SWE_LIBS) -o $(EXEC_DIR)/$(EXEC)

.PHONY: clean 
clean:
	@$(MAKE) -C $(MAIN_DIR) clean
	@$(MAKE) -C $(IO_DIR) clean
	@$(MAKE) -C $(RCM_DIR) clean
	@$(MAKE) -C $(LAYOUT_DIR) clean
	@$(MAKE) -C $(MPI_DIR) clean
	@$(MAKE) -C $(OCL_DIR) clean
	rm -f $(SOURCE_DIR)/*/*.a $(EXEC_DIR)/$(EXEC)

