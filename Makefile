OPENGL_LIB=-lGL -lGLU -lglut
EIGEN_LIB=-I /usr/include/eigen3

all: main

main: main.o particle.o mpm.o grid.o
	g++ -o main main.o particle.o mpm.o grid.o $(OPENGL_LIB) $(EIGEN_LIB)

main.o: main.cpp
	g++ -c main.cpp $(OPENGL_LIB) $(EIGEN_LIB)

mpm.o: mpm.cpp mpm.hpp 
	g++ -c mpm.cpp $(OPENGL_LIB) $(EIGEN_LIB)

particle.o: particle.cpp particle.hpp
	g++ -c particle.cpp $(OPENGL_LIB) $(EIGEN_LIB)

grid.o: grid.cpp grid.hpp
	g++ -c grid.cpp $(OPENGL_LIB) $(EIGEN_LIB)

.PHONY: clean
clean:
	rm -rf *.o test