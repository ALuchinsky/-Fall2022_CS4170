#CC = g++
UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
	CC = g++-11
else
	CC = g++
endif
LIBS = 
CFLAGS = -O2 -g -Wall -c -fmessage-length=0 
OPENMP = -fopenmp

RM := rm -rf

TARGET  = "OpenMP"

INCLUDE = -I../src

SUBDIRS  	 := $(wildcard ../) $(wildcard ../*/)
CPP_SRCS 	 := $(wildcard ../*/*.cpp) $(wildcard ../*/*/*.cpp) 

OBJS     := $(patsubst ../%.cpp, ./%.o, $(CPP_SRCS))
CPP_DEPS := $(patsubst ../%.cpp, ./%.d, $(CPP_SRCS))


# All Target
all: serial_Jacoby.exe parallel_Jacoby.exe

serial_Jacoby.exe: ../src/serial_Jacoby.cpp
	$(CC) -o $@ $< ../src/CStopWatch.cpp -I ../src/ $(OPENMP) $(ARCH) $(LIBS)

parallel_Jacoby.exe: ../src/parallel_Jacoby.cpp
	$(CC) -o $@ $< ../src/CStopWatch.cpp -I ../src/ $(OPENMP) $(ARCH) $(LIBS)


# Other Targets
clean:
	-$(RM) $(OBJS)$(CPP_DEPS)$(CPP_DEPS) $(TARGET)
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY: