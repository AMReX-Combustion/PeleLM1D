CANTERA_HOME ?= /home/marc/lib/cantera
#
# Note:  . ${CANTERA_HOME}/bin/setup_cantera to get the defs reqd below
#
include ${CANTERA_HOME}/include/cantera/Cantera.mak

CC=gcc
CXX=g++
RM=rm -f
CXXFLAGS=-pthread -O3 -Wno-inline -g -std=c++0x
CPPFLAGS=$(CANTERA_INCLUDES)
LDFLAGS=
LDLIBS=$(CANTERA_LIBS)

#EXEC=rpd
#EXEC=rpd0D
EXEC=rpd1Dopp

SRCS=$(EXEC).cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: $(EXEC)

clean:
	$(RM) $(OBJS) $(EXEC)

dist-clean: clean
	$(RM) *~ 
