#Nonlinear code makefile

all: run

run: testing
	./test

CC:=g++
CXX:=g++
ARMA:=-O2
CXXFLAGS:=-std=c++11 -fPIC -DANNELOGGER=nonllogger -I/usr/include -I/usr/include/python2.7 -Wall -Wno-unused-variable -Wno-return-local-addr $(ARMA)
LDFLAGS:=-Wextra -lpython2.7 -llapack -lblas -larmadillo -L/usr/lib -L/usr/lib/python2.7
RM:=rm -f

VPATH=%.cpp common common/bessel var

SOURCES=test.cpp globalconf.cpp tube.cpp material.cpp logger.cpp var.cpp gasvar.cpp
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
