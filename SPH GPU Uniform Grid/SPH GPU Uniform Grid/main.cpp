#include <glew.h>
#include <wglew.h>
#include <freeglut.h>

#include "sph.h"

const int number_of_particles = 8192 * 2;

const int dam_width = 256;     

const int scale = 2;    // Manipulates environment size, cell size, zoom of camera and point size.

SPH* sph;

// Called when window needs to be redrawn.
void display()
{
    // Clear screen.
    glClear(GL_COLOR_BUFFER_BIT);

    // Render particles.
    glDrawArrays(GL_POINTS, 0, number_of_particles);

    // Double buffering.
    glutSwapBuffers();
}

// Continuously called when events are not being received.
void idle()
{
    // Update function.
    sph->update();

    // Schedules a call to display function once all events are processed.
    glutPostRedisplay();
}

// Handles user input.
void keyboard(unsigned char key, int mouse_x, int mouse_y)
{
    // Setting direction of gravity in SPH simulation.
    switch (key)
    {
    case 'a':
        sph->setGravity(SPH::DIRECTION::LEFT);
        break;
    case 'd':
        sph->setGravity(SPH::DIRECTION::RIGHT);
        break;
    case 'w':
        sph->setGravity(SPH::DIRECTION::UP);
        break;
    case 's':
        sph->setGravity(SPH::DIRECTION::DOWN);
        break;
    }
}

int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(1024, 512);
    glutCreateWindow("GPU Uniform Grid. 8192 Particles.");

    glewInit();

    // Disable vertical sync
    wglSwapIntervalEXT(0);

    glMatrixMode(GL_PROJECTION);
    gluPerspective(60.0, 1024 / 512, 0.1f, 10000.0f);
    
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glPointSize(max(2, 4 / scale));
    glColor3f(0, 0.5, 1);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(-4.0f * scale, -1.0f * scale, -448.0f * scale);

    sph = new SPH(number_of_particles, dam_width, scale);

    // Defines an array of vertex data. 
    glVertexPointer(2, GL_FLOAT, 0, 0);

    // Enables use of vertex arrays. (needed to be able to call glDrawArrays)
    glEnableClientState(GL_VERTEX_ARRAY);

    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboard);

    glutMainLoop();

    delete sph;

    return 0;
}

