OPENGL_LIB=-lGL -lGLU -lglut
EIGEN_LIB=-I /usr/include/eigen3

all: main movie

main: main.o particle.o mpm.o grid.o material.o
	g++ -o main main.o particle.o mpm.o grid.o material.o $(OPENGL_LIB) $(EIGEN_LIB)

main.o: main.cpp
	g++ -c main.cpp $(OPENGL_LIB) $(EIGEN_LIB)

mpm.o: mpm.cpp mpm.hpp 
	g++ -c mpm.cpp $(OPENGL_LIB) $(EIGEN_LIB)

particle.o: particle.cpp particle.hpp
	g++ -c particle.cpp $(OPENGL_LIB) $(EIGEN_LIB)

grid.o: grid.cpp grid.hpp
	g++ -c grid.cpp $(OPENGL_LIB) $(EIGEN_LIB)

material.o: material.cpp material.hpp
	g++ -c material.cpp $(OPENGL_LIB) $(EIGEN_LIB)

movie: movie.o 
	g++ -o movie movie.o $(OPENGL_LIB) $(EIGEN_LIB)

movie.o: movie.cpp 
	g++ -c movie.cpp $(OPENGL_LIB) $(EIGEN_LIB)

.PHONY: clean
clean:
	rm -rf *.o test