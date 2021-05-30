#ifndef MDSYSTEMGL_H
#define MDSYSTEMGL_H

#include "MDSystem.h"
#include <QOpenGLFunctions>
#include <QtOpenGLExtensions/QOpenGLExtensions>

class MDSystemGL : protected QOpenGLFunctions, QOpenGLExtension_ARB_multitexture
{
public:
    bool m_renderSpheres;

    // Display lists for rendering the spheres
    int dListSphere;
    bool isDList;

    bool usepbo;

    // Rendering stuff
    unsigned int m_program, m_program2;
    unsigned int m_vertexShader, m_pixelShader, m_texture;

    int ox, oy;
    float camera_trans[3];
    float camera_rot[3];
    float camera_trans_lag[3];
    float camera_rot_lag[3];
    int buttonState;

    MDSystem *m_system;

    MDSystemGL(const MDSystem::MDSystemConfiguration& config = MDSystem::MDSystemConfiguration(), bool renderSpheres = true);
    ~MDSystemGL(void) { delete m_system; }

    void Reinitialize(const MDSystem::MDSystemConfiguration& config = MDSystem::MDSystemConfiguration(), bool renderSpheres = true);

    void draw();
    void InitGL();
    void printrdf();
    void _createTexture(int resolution);
    void initSphere();

    void setRenderMode(bool renderSpheres) { m_renderSpheres = renderSpheres; }
    void setHardwareMode(bool useCUDA) { m_system->setHardwareMode(useCUDA); }
    void setCanonical(bool canonical) { m_system->setCanonical(canonical); }
    void setPeriodicBoundaryCondition(bool periodic) { m_system->setPeriodicBoundaryCondition(periodic); }
};

#endif // MDSYSTEMGL_H
