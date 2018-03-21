CXXFLAGS=-std=c++14

INCLUDES = -I/home/rahulgrit/miniconda2/include/

LDLFLAGS = -L/home/rahulgrit/miniconda2/lib/

TARGETS = invtmat

all: $(TARGETS)

clean:
	rm -f $(TARGETS)

distclean: clean
	rm -f *.o *~

invtmat: invtmat.cpp matinv.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDLFLAGS) -o $@ $<

