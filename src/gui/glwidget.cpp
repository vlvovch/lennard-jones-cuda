#include <QtGui>
//#include <QtOpenGL>

#include <math.h>

#include "glwidget.h"
#define PI 3.14159265

#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE  0x809D
#endif

//! [0]
GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
    Initialize();
}

GLWidget::GLWidget(const QGLFormat & format, QWidget *parent)
    : QGLWidget(format, parent)
{
    Initialize();
}
//! [0]

//! [1]
GLWidget::~GLWidget()
{
}

void GLWidget::Initialize()
{
    timer.start();
    anim = new QTimer(this);
    connect(anim, SIGNAL(timeout()), this, SLOT(change()));
    //anim->start(10);
    wrk = true;
    fpsTimer.start();
    frames = 0;
    fps = 0.;
    setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    whint = 800;
    hhint = 600;
    xRot = yRot = zRot = 0;

    video_frame = 0;

    record_frames = false;
}
//! [1]

//! [2]
QSize GLWidget::minimumSizeHint() const
{
    return QSize(50, 50);
}
//! [2]

//! [3]
QSize GLWidget::sizeHint() const
//! [3] //! [4]
{
    return QSize(600, 370);
    //return QSize(whint, hhint);
}
//! [4]

static void qNormalizeAngle(int &angle)
{
    while (angle < 0)
        angle += 360 * 16;
    while (angle > 360 * 16)
        angle -= 360 * 16;
}

//! [5]
void GLWidget::setXRotation(int angle)
{
    //qNormalizeAngle(angle);
    if (angle != xRot) {
        xRot = angle;
        systGL->camera_rot[0] = angle / 5.;
        emit xRotationChanged(angle);
        //updateGL();
    }
}
//! [5]

void GLWidget::setYRotation(int angle)
{
    //qNormalizeAngle(angle);
    if (angle != yRot) {
        yRot = angle;
        systGL->camera_rot[1] = angle / 5.;
        emit yRotationChanged(angle);
        //updateGL();
    }
}

void GLWidget::setZRotation(int angle)
{
    //qNormalizeAngle(angle);
    if (angle != zRot) {
        zRot = angle;
        systGL->camera_rot[2] = angle / 5.;
        emit zRotationChanged(angle);
        //updateGL();
    }
}

void GLWidget::change()
{
    if (!wrk) return;
    for(int i=0;i<pre;++i) 
      systGL->m_system->Integrate(dt);
      //systGL->m_system->Integrate(dt/pre);
    updateGL();
    emit updateGraph();

    
    if (record_frames)
      this->grabFrameBuffer().save("out_anim/screen" + QStringLiteral("%1").arg(video_frame++, 5, 10, QLatin1Char('0')) + ".png");
	  //printf("%lf\n", (systGL->m_system->K + systGL->m_system->V)/ systGL->m_system->m_config.N);

}

void GLWidget::changemode(int mode)
{
    //syst->changedrawmode((bool)mode);
}

//! [6]
void GLWidget::initializeGL()
{
    systGL->InitGL();
}
//! [6]

//! [7]
void GLWidget::paintGL()
{
    //glMatrixMode(GL_MODELVIEW);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //glLoadIdentity();

    //glMatrixMode(GL_PROJECTION);
    //glLoadIdentity();
    //perspectiveGL(60.0, (double)(width())/height(), 1.0, 600.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if (wrk) systGL->draw();
    glColor3ub(255,255,255);
    glColor3ub(0,0,0);
    frames++;
    if (frames%20==0)
    {
        fps = frames * 1.e3 / fpsTimer.restart();
        frames = 0;
    }
    double u = (systGL->m_system->K + systGL->m_system->V) / systGL->m_system->m_config.N;
    double u_av = systGL->m_system->av_U_tot / systGL->m_system->av_iters / systGL->m_system->m_config.N;
    double T = systGL->m_system->T;
    double T_av = systGL->m_system->av_T_tot / systGL->m_system->av_iters;
    double Z    = systGL->m_system->P / T_av / systGL->m_system->m_config.rho;
    double Z_av = (systGL->m_system->av_p_tot / systGL->m_system->av_iters) / T_av / systGL->m_system->m_config.rho;

    //return;

    //renderText(10, 20, QString(tr("%1").arg(syst->x[0])));
    //renderText(10, height() - (60 * height())/500, QString(tr("Energy per particle: %1 meV")).arg((systGL->m_system->K + systGL->m_system->V) * ent / eneV * 1.e3 / systGL->m_system->m_config.N));
    //renderText(10, height() - (40 * height())/500, QString(tr("Instant temperature: %1 K")).arg((systGL->m_system->T) * ent / kb));
    renderText(10, height() - (60 * height()) / 500, QString(tr("Energy per particle: u* = %1, <u*> = %2")).arg((systGL->m_system->K + systGL->m_system->V) / systGL->m_system->m_config.N).arg(u_av));
    renderText(10, height() - (40 * height()) / 500, QString(tr("Instant temperature: T* = %1, <T*> = %2")).arg((systGL->m_system->T)).arg(T_av));
    //renderText(10, height() - (20 * height()) / 500, QString(tr("Pressure: p* = %1, <p*> = %2")).arg((systGL->m_system->P)));
    renderText(10, height() - (20 * height()) / 500, QString(tr("Compressibility: Z* = %1, <Z*> = %2")).arg(Z).arg(Z_av));
    renderText(10, (20 * height())/500, QString(tr("FPS: %1")).arg(fps));
    renderText(10, (40 * height()) / 500, QString(tr("Time: %1")).arg(systGL->m_system->getTime()));
    //renderText(10, (40 * height())/500, QString(tr("Time: %1 ps")).arg(systGL->m_system->getTime() * timt * 1.e12));

    //glDisable(GL_MULTISAMPLE);
    /*QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    p->drawtext(&painter, width(), height());
    painter.end();*/
    //swapBuffers();
}
//! [7]

//! [8]
void GLWidget::resizeGL(int width, int height)
{
    if (height==0)										// Prevent A Divide By Zero By
    {
        height=1;										// Making Height Equal One
    }

    glViewport(0,0,width,height);						// Reset The Current Viewport

    glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
    glLoadIdentity();									// Reset The Projection Matrix

    // Calculate The Aspect Ratio Of The Window
    /*gluPerspective*/perspectiveGL(60.0f,(GLfloat)width/(GLfloat)height,0.1f,500.0f);

    glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
    glLoadIdentity();	            								// Reset The Modelview Matrix
}
//! [8]

//! [9]
void GLWidget::mousePressEvent(QMouseEvent *event)
{
    lastPos = event->pos();
}
//! [9]

//! [10]
void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    int dx = event->x() - lastPos.x();
    int dy = event->y() - lastPos.y();

    if (event->buttons() & Qt::LeftButton) {
        setXRotation(xRot + dy);
        setYRotation(yRot + dx);
    } else if (event->buttons() & Qt::RightButton) {
        /*setXRotation(xRot + dy);
        setZRotation(zRot + dx);*/
        systGL->camera_trans[2] += (dy / 10.) * 5;
    }
    lastPos = event->pos();
    updateGL();
}
//! [10]
