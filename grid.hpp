#pragma once

#include "Eigen/Dense"
#include "GL/gl.h"
#include "GL/glut.h"

class Grid
{
public:
    void render();

    Eigen::Vector2d position;

    // Intermediate
    double mass;
    Eigen::Vector2d linearMomentum;
    Eigen::Vector2d newLinearMomentum;

    Eigen::Vector2d force;
    Eigen::Vector2d newPosition;
};