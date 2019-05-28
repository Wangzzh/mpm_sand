#include <cstdlib>
#include <iostream>

#include "GL/gl.h"
#include "GL/glut.h"

#include "particle.hpp"
#include "grid.hpp"
#include "mpm.hpp"

MPM* mpm;

bool simulate = false;

int n = 30;
double dTime = 0.004;

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
}

void idle() {
    if (simulate) {
        std::cout << "Auto step" << std::endl;
        mpm->step();
    }
    glutPostRedisplay();
}

int main(int argc, char* argv[]) {
    MaterialParameters material = MaterialParameters(140000, 0.2, 10, 0.025, 0.0075);
    mpm = new MPM(n, dTime, material);

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