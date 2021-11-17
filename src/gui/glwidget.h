#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include <QElapsedTimer>
#include <QTimer>
#include "constants.h"
#include "MDSystemGL.h"
//! [0]
class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget(QWidget *parent = 0);
    GLWidget(const QGLFormat & format, QWidget *parent = 0);
    ~GLWidget();

    void Initialize();

    QSize minimumSizeHint() const;
    QSize sizeHint() const;

    MDSystemGL *systGL;

    QTimer *anim;

    bool wrk;

    // Time step between two frames
    double dt;
    // Numerical integration per time step
    int pre;

    bool record_frames;


    int whint, hhint;
//! [0]

//! [1]
public slots:
    void setXRotation(int angle);
    void setYRotation(int angle);
    void setZRotation(int angle);
    void change();
    void changemode(int);

signals:
    void xRotationChanged(int angle);
    void yRotationChanged(int angle);
    void zRotationChanged(int angle);
    void updateGraph();
//! [1]

//! [2]
protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
//! [2]

//! [3]
private:
    int xRot;
    int yRot;
    int zRot;
    QPoint lastPos;
    QElapsedTimer timer;
    QElapsedTimer fpsTimer;
    int frames;
    double fps;

    int video_frame;
};
//! [3]

#endif // GLWIDGET_H
