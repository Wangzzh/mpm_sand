#include "particle.hpp"

bool Particle::renderV = false;
bool Particle::renderF = false;
bool Particle::renderFE = false;
bool Particle::renderB = false;

void Particle::render() {
    glColor3f(color(0), color(1), color(2));
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glRectd(position(0) - 0.003, position(1) - 0.003, position(0) + 0.003, position(1) + 0.003);
    // glTranslatef(position[0], position[1], 0.);
    // glutSolidSphere(0.004, 3, 3);
    glPopMatrix();

    if (Particle::renderV) {
        glColor3f(1, 1, 0);
        double velocityFactor = 0.02;
        glBegin(GL_LINES);
        glVertex3f(position[0], position[1], 0.0);
        glVertex3f(position[0] + velocity[0] * velocityFactor, 
            position[1] + velocity[1] * velocityFactor, 0.);
        glEnd();
    }
    
    if (Particle::renderFE) {
        double velocityFactor = 0.05;
        glBegin(GL_LINES);
        glColor3f(0.5, 1, 1);
        glVertex3f(position[0], position[1], 0.0);
        glVertex3f(position[0] + FE(0,0) * velocityFactor, 
            position[1] + FE(1,0) * velocityFactor, 0.);
        glColor3f(1, 0.5, 1);
        glVertex3f(position[0], position[1], 0.0);
        glVertex3f(position[0] + FE(1,0) * velocityFactor, 
            position[1] + FE(1,1) * velocityFactor, 0.);
        glEnd();
    }
    
    if (Particle::renderB) {
        double velocityFactor = 0.05;
        glBegin(GL_LINES);
        glColor3f(0.5, 0, 0.5);
        glVertex3f(position[0], position[1], 0.0);
        glVertex3f(position[0] + B(0,0) * velocityFactor, 
            position[1] + B(1,0) * velocityFactor, 0.);
        glColor3f(0., 0.5, 0.5);
        glVertex3f(position[0], position[1], 0.0);
        glVertex3f(position[0] + B(1,0) * velocityFactor, 
            position[1] + B(1,1) * velocityFactor, 0.);
        glEnd();
    }

    if (Particle::renderF) {
        double velocityFactor = 0.05;
        Eigen::Matrix2d F = FE * FP;
        glBegin(GL_LINES);
        glColor3f(0.8, 1, 1);
        glVertex3f(position[0], position[1], 0.0);
        glVertex3f(position[0] + F(0,0) * velocityFactor, 
            position[1] + F(1,0) * velocityFactor, 0.);
        glColor3f(1, 0.8, 1);
        glVertex3f(position[0], position[1], 0.0);
        glVertex3f(position[0] + F(1,0) * velocityFactor, 
            position[1] + F(1,1) * velocityFactor, 0.);
        glEnd();
    }
}

void Particle::calculateWeights() {
    xWeight.clear(); yWeight.clear();
    xWeightGradient.clear(); yWeightGradient.clear();

    double x = xDiff;
    xWeight.push_back(-(1+x) * (1+x) * (1+x) / 6. + (1+x) * (1+x) - 2 * (1+x) + 4./3.);
    xWeight.push_back(0.5 * x * x * x - x * x + 2./3.);
    xWeight.push_back(0.5 * (1-x) * (1-x) * (1-x) - (1-x) * (1-x) + 2./3.);
    xWeight.push_back(-(2-x) * (2-x) * (2-x) / 6. + (2-x) * (2-x) - 2 * (2-x) + 4./3.);
    
    xWeightGradient.push_back(-(1+x) * (1+x) / 2. + 2 * (1+x) - 2);
    xWeightGradient.push_back(1.5 * x * x - 2 * x);
    xWeightGradient.push_back(-1.5 * (1-x) * (1-x) + 2 * (1-x));
    xWeightGradient.push_back((2-x) * (2-x) / 2. - 2 * (2-x) + 2);
    
    x = yDiff;
    yWeight.push_back(-(1+x) * (1+x) * (1+x) / 6. + (1+x) * (1+x) - 2 * (1+x) + 4./3.);
    yWeight.push_back(0.5 * x * x * x - x * x + 2./3.);
    yWeight.push_back(0.5 * (1-x) * (1-x) * (1-x) - (1-x) * (1-x) + 2./3.);
    yWeight.push_back(-(2-x) * (2-x) * (2-x) / 6. + (2-x) * (2-x) - 2 * (2-x) + 4./3.);
    
    yWeightGradient.push_back(-(1+x) * (1+x) / 2. + 2 * (1+x) - 2);
    yWeightGradient.push_back(1.5 * x * x - 2 * x);
    yWeightGradient.push_back(-1.5 * (1-x) * (1-x) + 2 * (1-x));
    yWeightGradient.push_back((2-x) * (2-x) / 2. - 2 * (2-x) + 2);

}