CXX=clang++
CXXFLAGS=-std=c++11 -Iinclude/ -I/usr/lib/gcc/x86_64-pc-linux-gnu/10.2.0/include/ -O3
LDFLAGS=lib/libpng16.a lib/libz.a -no-pie -lpthread 

OBJECTS=src/main.o src/image.o src/stencil.o

stencil: $(OBJECTS)
	$(CXX) -o bin/stencil $(OBJECTS) $(LDFLAGS)

all:	stencil

run:	all
	bin/stencil test-image.png

clean:
	rm -f $(OBJECTS) bin/stencil data.txt src/*~ *~
