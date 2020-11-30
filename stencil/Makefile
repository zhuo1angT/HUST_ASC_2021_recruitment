CXX=clang++
CXXFLAGS=-std=c++11 -Iinclude/ -I/usr/lib/gcc/x86_64-pc-linux-gnu/10.2.0/include/ -O3
LDFLAGS=lib/libpng16.a lib/libz.a -no-pie -lpthread 

BASELINEOBJS=src/main.o src/image.o src/stencil.o
MYOBJS=src/main.o src/image.o src/mine.o src/stencil.o

clean:
	rm -f $(OBJECTS) bin/stencil data.txt src/*~ *~
	rm -f $(MYOBJS) 

org: $(BASELINEOBJS)
	$(CXX) -o bin/stencil $(BASELINEOBJS) $(LDFLAGS)

mine: $(MYOBJS)
	$(CXX) -o bin/stencil $(MYOBJS) $(LDFLAGS)

check_org: org
	bin/stencil IPCC.png

check_mine: mine
	bin/stencil IPCC.png