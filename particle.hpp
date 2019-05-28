#pragma once

#include <iostream>
#include <vector>

#include "Eigen/Dense"
#include "GL/gl.h"
#include "GL/glut.h"

class Particle
{
public:
    void render();

    double mass;
    double density;
    double volume;

    Eigen::Vector2d position;
    Eigen::Vector2d velocity;
    
    Eigen::Matrix2d FE = Eigen::Matrix2d::Identity();
    Eigen::Matrix2d FP = Eigen::Matrix2d::Identity();

    // intermediate results
    int xLeft, yLeft;
    double xDiff, yDiff;
    std::vector<double> xWeight, yWeight;
    std::vector<double> xWeightGradient, yWeightGradient;

    void calculateWeights();
};