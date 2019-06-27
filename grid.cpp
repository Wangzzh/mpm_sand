#include "grid.hpp"

bool Grid::renderV = false;
bool Grid::renderF = false;

void Grid::render() {
    glColor3f(color(0), color(1), color(2));
    glMatrixMode(GL_MODELVIEW);
    // glPushMatrix();
    // glLoadIdentity();
    glRectd(position(0) - 0.002, position(1) - 0.002, position(0) + 0.002, position(1) + 0.002);
    // glTranslatef(position[0], position[1], 0.);
    // glutSolidSphere(0.002, 3, 3);
    // glPopMatrix();
    
    if (Grid::renderV) {
        glColor3f(1, 0, 0);
        double velocityFactor = 2;
        glBegin(GL_LINES);
        glVertex3f(position[0], position[1], 0.0);
        glVertex3f(position[0] + linearMomentum[0] * velocityFactor, 
            position[1] + linearMomentum[1] * velocityFactor, 0.);
        glEnd();
    }
    
    if (Grid::renderF) {
        glColor3f(0, 1, 0);
        double velocityFactor = 0.5;
        glBegin(GL_LINES);
        glVertex3f(position[0], position[1], 0.0);
        glVertex3f(position[0] + force[0] * velocityFactor, 
            position[1] + force[1] * velocityFactor, 0.);
        glEnd();
    }
}