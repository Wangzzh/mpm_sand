#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>

#include "GL/gl.h"
#include "GL/glut.h"

bool play = false;

std::ifstream f;
std::string filename = "out.txt";

double dTime;

void playFrame() {
    std::string line;
    if (std::getline(f, line)) {
        std::istringstream iss(line);
        double x, y;
        while (iss >> x >> y) {
            double r,g,b;
            iss >> r >> g >> b;
            glColor3f(r, g, b);
            glMatrixMode(GL_MODELVIEW);
            glPushMatrix();
            glLoadIdentity();
            glRectd(x - 0.003, y - 0.003, x + 0.003, y + 0.003);
            glPopMatrix();
        }
    } else {
        std::string line;
        f.close();
        f.open(filename);
        std::getline(f, line);
    }
    glFlush();
}

void display() {
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);

    if (play) {
        playFrame();
    }

}

void keyboard(unsigned char key, int x, int y) {
    if (key == ' ') {
        play = !play;
        std::cout << "Play: " << play << std::endl;
    }
    if (key == 'g' && !play) {
        playFrame();
    }
}

void idle() {
    glutPostRedisplay();
}

int main(int argc, char* argv[]) {
    if (argc >= 2) {
        filename = std::string(argv[1]);
    }

    f.open(filename);
    std::string line;
    std::getline(f, line);
    std::cout << "Dtime " << line << std::endl;
    std::istringstream iss(line);
    iss >> dTime;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("MPM Movie");
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 1, 0, 1, 0, 100);

    glutMainLoop();
    return 0;
}
