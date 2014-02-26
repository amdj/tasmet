#Nonlinear code makefile

all: run

run: testing
	./test

CC:=g++
CXX:=g++
ARMA:=-O2
CXXFLAGS:=-std=c++11 -fPIC -DANNELOGGER=nonllogger -I/usr/include -I/usr/include/python2.7 -Wall $(ARMA)
LDFLAGS:=-Wextra -lpython2.7 -larmadillo -L/usr/lib -L/usr/lib/python2.7
RM:=rm -f

VPATH=%.cpp common common/bessel

SOURCES=test.cpp globalconf.cpp material.cpp logger.cpp
HEADERS: vtypes.h
OBJS:=$(addprefix .obj/,$(notdir $(SOURCES:.cpp=.o)))




testing: $(OBJS)
	$(CXX) $(CXXFLAGS) -o test $(OBJS) $(LDFLAGS)

$(OBJS): | .obj

.obj/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


$(SOURCES): $(HEADERS)

.obj:
	mkdir .obj

.PHONY: clean
clean:
	$(RM) -rf .obj test
