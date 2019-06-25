#include <cstdlib>
#include <iostream>

#include "GL/gl.h"
#include "GL/glut.h"

#include "particle.hpp"
#include "grid.hpp"
#include "mpm.hpp"

MPM* mpm;

bool simulate = false;
bool pour = false;
int pour_interval = 4;
int pour_counter = 0;

int n = 100;
double dTime = 0.0005;

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
    if (key == 'p') {
        pour = !pour;
        std::cout << "Pouring " << pour << std::endl;
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
        if (pour) {
            pour_counter ++;
            if (pour_counter % pour_interval == 0) {
                mpm -> addCube(Eigen::Vector2d(0.5, 0.3), Eigen::Vector2d(0.05, 0.01), 0, 1, 1, mpm->materials[0], Eigen::Vector3d(1, 1, 1));
                pour_counter -= pour_interval;
            }
        }
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