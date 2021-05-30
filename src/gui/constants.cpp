#include "constants.h"
#include <QGLWidget>
#include <cmath>

double kb = 1.3806488e-23;
double mu = 1.660538921e-27;
double eneV = 1.602176565e-19;

double rt = 0.34e-9;
double ent = 120*kb;
double masst = 39.948; // Argon
double timt = rt/sqrt(ent/(mu*masst));

// Replaces gluPerspective. Sets the frustum to perspective mode.
// fovY     - Field of vision in degrees in the y direction
// aspect   - Aspect ratio of the viewport
// zNear    - The near clipping distance
// zFar     - The far clipping distance

void perspectiveGL( GLdouble fovY, GLdouble aspect, GLdouble zNear, GLdouble zFar )
{
    const GLdouble pi = 3.1415926535897932384626433832795;
    GLdouble fW, fH;
    fH = tan( fovY / 360 * pi ) * zNear;
    fW = fH * aspect;
    glFrustum( -fW, fW, -fH, fH, zNear, zFar );
}


