
# This is the default mode to run gNB code,
# When no command line flag is specified
MODE = RUN_GNB_PHY

csrc =  $(wildcard *.c)

CC = gcc

#obj = $(csrc:.c=.o)

# Command line inputs while compiling
# Example: make -B MODE=PROFILE_GNB

ifeq "$(MODE)" "DEBUG"
	MODE_FLAGS = -DDEBUG

else ifeq "$(MODE)" "PROFILE_GNB"
	MODE_FLAGS = -DPROFILE_GNB

else ifeq "$(MODE)" "RUN_GNB_PHY"
	MODE_FLAGS = -DRUN_GNB_PHY

endif

all : clean wnNrPhyExe

wnNrPhyExe: #$(obj)
	@echo "gNB PHY is running with MODE : $(MODE)"
	@echo "Compiling with : $(CC)"
	$(CC) $(MODE_FLAGS) $(CFLAGS) -o wnNrPhyExe $(csrc) $(CFLAGS) $(LDFLAGS)
	#$(CC) -o $@ $^ $(CFLAGS) -o wnNrPhyExe

#-Werror
CFLAGS += -DAVX_2 -O1
CFLAGS += -g -lm -unroll-aggressive -Winline -D_GNU_SOURCE -mavx2 -mfma -lpthread  -lfftw3f -std=c99 -O3 
CFLAGS += $(WERROR_FLAGS)
LDFLAGS  = -lm

.PHONY: clean
clean:
	rm -rf $(obj) wnNrPhyExe *.txt
	@echo "All object files, text files and executable have been removed"
