#ifndef CONSTANTS_H
#define CONSTANTS_H

// Physical constants

// Boltzmann constant
extern double kb;

// Atomic mass unit
extern double mu;

// 1 electronvolt in joules
extern double eneV;


// Argon parameters

// Hard-core radius
extern double rt;

// Attractive well
extern double ent;

// Atomic mass
extern double masst;

// Conversion of time unit to seconds
extern double timt;


// Replaces gluPerspective. Sets the frustum to perspective mode.
// fovY     - Field of vision in degrees in the y direction
// aspect   - Aspect ratio of the viewport
// zNear    - The near clipping distance
// zFar     - The far clipping distance
void perspectiveGL( double fovY, double aspect, double zNear, double zFar );


#endif // CONSTANTS_H
