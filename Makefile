#Nonlinear code makefile

all: run

run: testing
	./test
test_drag: test_dragops
	./test_drag

CC:=g++
CXX:=g++
ARMA:=-O2
CXXFLAGS:=-std=c++11 -fPIC -DANNELOGGER=nonllogger -I/usr/include -I/usr/include/python2.7 -I/usr/local/include -Wall -Wno-unused-variable -Wno-unused-but-set-variable -Wno-return-local-addr $(ARMA)
CXXFLAGS+=-MD -MP
LDFLAGS:= -L/usr/local/lib -L/usr/lib -L/usr/lib/python2.7 -Wextra -lpython2.7 -lmath_common -llapack -lblas -larmadillo
RM:=rm -f
MAKEDEPEND=gcc -M $(CXXFLAGS)

VPATH=%.cpp var tube

SOURCES=globalconf.cpp tube.cpp vertex.cpp var.cpp gasvar.cpp geom.cpp drag.cpp  continuityeq.cpp momentumeq.cpp tubeequation.cpp energyeq.cpp stateeq.cpp solidenergyeq.cpp
HEADERS: vtypes.h
OBJS:=$(addprefix .obj/,$(notdir $(SOURCES:.cpp=.o)))


test_dragops: $(OBJS) .obj/test_drag.o
	$(CXX) $(CXXFLAGS) -o test_drag .obj/test_drag.o $(OBJS) $(LDFLAGS)

testing: $(OBJS) .obj/test.o
	$(CXX) $(CXXFLAGS) -o test .obj/test.o $(OBJS) $(LDFLAGS)

$(OBJS): | .obj


#INCLUDING AUTO-GENERATION OF DEPENDENCY FILES
.obj/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MD -c $< -o $@
	@cp .obj/$*.d .obj/$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
          -e '/^$$/ d' -e 's/$$/ :/' < .obj/$*.d >> .obj/$*.P; \
           rm -f .obj/$*.d


$(SOURCES): $(HEADERS)

.obj:
	mkdir .obj

.PHONY: clean
clean:
	$(RM) -rf .obj test test_drag

-include $(OBJS:.o=.P)


