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

TOOLS_FILES = ../src/CStopWatch.cpp ../src/utils.cpp

# All Target
all: serial_Jacoby.exe parallel_Jacoby.exe Power.exe Household.exe

Household.exe: ../src/Household.cpp
	$(CC) -o $@ $< ${TOOLS_FILES} -I ../src/ $(OPENMP) $(ARCH) $(LIBS)


serial_Jacoby.exe: ../src/serial_Jacoby.cpp
	$(CC) -o $@ $< ${TOOLS_FILES} -I ../src/ $(OPENMP) $(ARCH) $(LIBS)

parallel_Jacoby.exe: ../src/parallel_Jacoby.cpp
	$(CC) -o $@ $< ${TOOLS_FILES} -I ../src/ $(OPENMP) $(ARCH) $(LIBS)

Power.exe: ../src/Power.cpp
	$(CC) -o $@ $< ${TOOLS_FILES} -I ../src/ $(OPENMP) $(ARCH) $(LIBS)


# Other Targets
clean:
	rm *.exe

.PHONY: all clean dependents
.SECONDARY: