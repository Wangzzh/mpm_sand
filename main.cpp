#include <cstdlib>
#include <iostream>

#include "GL/gl.h"
#include "GL/glut.h"

#include "particle.hpp"
#include "grid.hpp"
#include "mpm.hpp"

MPM* mpm;

bool simulate = false;

int n = 50;
double dTime = 0.002;

void display() {
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);

    mpm->render();
    // std::cout << "step" << std::endl;

    glFlush();
}

void keyboard(unsigned char key, int x, int y) {
    if (key == 'g'  && !simulate) {
        std::cout << "Step" << std::endl;
        mpm -> step();
        // mpm -> render();
    }
    if (key == ' ') {
        simulate = !simulate;
    }
    if (key == '1') {
        Particle::renderV = !Particle::renderV;
    }
    if (key == '2') {
        Particle::renderFE = !Particle::renderFE;
    }
    if (key == '3') {
        Particle::renderF = !Particle::renderF;
    }
    if (key == '4') {
        Particle::renderB = !Particle::renderB;
    }
    if (key == '6') {
        Grid::renderV = !Grid::renderV;
    }
    if (key == '7') {
        Grid::renderF = !Grid::renderF;
    }
}

void idle() {
    if (simulate) {
        // std::cout << "Auto step" << std::endl;
        mpm->step();
    }
    glutPostRedisplay();
}

int main(int argc, char* argv[]) {
    mpm = new MPM(n, dTime);    

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("MPM Snow Simulation");
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 1, 0, 1, 0, 100);

    glutMainLoop();
    return 0;
}