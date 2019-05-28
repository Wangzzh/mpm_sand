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

class MaterialParameters 
{
public:
    MaterialParameters() {}
    MaterialParameters(double youngsModulus, double poissonsRatio, 
        double hardening, double criticalCompression, double criticalStretch);
    double E, nu, xsi, thetaC, thetaS, lambda, mu;
};

class MPM
{
public:
    MPM(int nGrid, double timeStep, MaterialParameters material);
    ~MPM();
    void render();
    void step();

private:
    std::vector<double> calculateWeights(double distToLeft);
    
    void particleToGrid();
    void computeParticleDensity();
    void computeGridForce();
    void computeGridVelocity();
    void updateDeformation();
    void updateParticleVelocity();
    void handleParticleCollision();
    void updateParticlePosition();

    double timeStep;
    double time;

    int nGrid;
    std::vector<Particle*> particles;
    std::vector<std::vector<Grid*>> grids;
    
    MaterialParameters material;
};