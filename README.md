# Material Point Simulation of various materials

C++ implementation of Material Point Method (MPM) Simulation of fluid-like and sand-like material. The implementation is highly based on "Drucker-prager elastoplasticity for sand animation" by Klár G et. al.

The latest version is on "sand" branch. 

### Dependencies

- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [OpenGL](https://www.opengl.org/)

### How to run

- Install the dependencies. 
- `make` to build the project. You may need to change the eigen library path in the Makefile.
- `./main` to run the simulation.

- If set `output=true` in main.cpp, a data file will be recorded. Run `./movie <data file name>` will produce a replay of the simualtion.

### Demonstration

A demo video is avaiable [HERE](https://youtu.be/Kck63ryrxHc).

### Reference
- Stomakhin A, Schroeder C, Chai L, et al. A material point method for snow simulation[J]. ACM Transactions on Graphics (TOG), 2013, 32(4): 102.
- Klár G, Gast T, Pradhana A, et al. Drucker-prager elastoplasticity for sand animation[J]. ACM Transactions on Graphics (TOG), 2016, 35(4): 103.
- Tampubolon A P, Gast T, Klár G, et al. Multi-species simulation of porous sand and water mixtures[J]. ACM Transactions on Graphics (TOG), 2017, 36(4): 105.
