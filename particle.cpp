#include "particle.hpp"

void Particle::render() {
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(position[0], position[1], 0.);
    glutSolidSphere(0.002, 10, 10);
    glPopMatrix();

    // double velocityFactor = 0.1;
    // glBegin(GL_LINES);
    // glVertex3f(position[0], position[1], 0.0);
    // glVertex3f(position[0] + velocity[0] * velocityFactor, 
    //     position[1] + velocity[1] * velocityFactor, 0.);
    // glEnd();
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
    xWeightGradient.push_back(1.5 * x * x - 2 * x + 2./3.);
    xWeightGradient.push_back(-1.5 * (1-x) * (1-x) + 2 * (1-x) + 2./3.);
    xWeightGradient.push_back((2-x) * (2-x) / 3. - 2 * (2-x) + 2);
    
    x = yDiff;
    yWeight.push_back(-(1+x) * (1+x) * (1+x) / 6. + (1+x) * (1+x) - 2 * (1+x) + 4./3.);
    yWeight.push_back(0.5 * x * x * x - x * x + 2./3.);
    yWeight.push_back(0.5 * (1-x) * (1-x) * (1-x) - (1-x) * (1-x) + 2./3.);
    yWeight.push_back(-(2-x) * (2-x) * (2-x) / 6. + (2-x) * (2-x) - 2 * (2-x) + 4./3.);
    
    yWeightGradient.push_back(-(1+x) * (1+x) / 2. + 2 * (1+x) - 2);
    yWeightGradient.push_back(1.5 * x * x - 2 * x + 2./3.);
    yWeightGradient.push_back(-1.5 * (1-x) * (1-x) + 2 * (1-x) + 2./3.);
    yWeightGradient.push_back((2-x) * (2-x) / 3. - 2 * (2-x) + 2);
}