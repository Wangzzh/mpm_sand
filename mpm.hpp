#pragma once

#include <cstdlib>
#include <cmath>
#include <ctime>
#include <random>
#include <vector>

#include "Eigen/Dense"
#include "Eigen/SVD"
#include "GL/gl.h"
#include "GL/glut.h"

#include "particle.hpp"
#include "grid.hpp"
#include "material.hpp"

class MPM
{
public:
    MPM(int nGrid, double timeStep);
    ~MPM();
    void render();
    void step();

    void particleToGrid();
    void computeParticleDensity();
    void computeGridForce();
    void computeGridVelocity();
    void updateDeformation();
    void updateParticleVelocity();
    void handleParticleCollision();
    void updateParticlePosition();

    void addCube(Eigen::Vector2d position, Eigen::Vector2d size, double angle, 
        int div, int random, MaterialParameters* material, Eigen::Vector3d color);

    double timeStep;
    double time;

    int nGrid;
    std::vector<Particle*> particles;
    std::vector<std::vector<Grid*>> grids;
    std::vector<MaterialParameters*> materials;
};