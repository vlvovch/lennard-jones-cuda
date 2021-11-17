#ifndef MAINWINDOW_H
#define MAINWINDOW_H

//#include <QtGui/QMainWindow>
#include <QtGui>
#include "glwidget.h"
#include "qcustomplot.h"


class MainWindow : public QMainWindow
{
    Q_OBJECT

protected:
    void keyPressEvent(QKeyEvent *event);
    void drawRDF();
    void drawvel();

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

    MDSystem::MDSystemConfiguration getConfig();

private:
    GLWidget *glWidget;

    QSpinBox *spinN, *spinpre;
    QDoubleSpinBox *spinrho, *spinT, *spindtau;

    QRadioButton *RBMicro, *RBCanon;
    QRadioButton *RBPeriodic, *RBWalls;
    QRadioButton *RBVerlet;
    QRadioButton *RBCPU, *RBCUDA;
    QRadioButton *RBSpheres, *RBPoints;
    QPushButton *buttonRestart, *buttonCont, *buttonStop, *buttonReaverage;

    QCheckBox* CBframes;

    QVBoxLayout *lhsLayout;

    QGridLayout *plotsLayout;
    QWidget *plotsWidget;

    QCustomPlot *ploten, *plotrdf, *plotvel;
    SplineFunction rdftot;
    int rdfiter;
private slots:
    void restart();
    void reav();
    void cont();
    void stop();
    void resizeEvent(QResizeEvent *event);
    void updateEner();
    void updateRenderMode();
    void updateEnsemble();
    void updateBoundary();
    void updateHardwareMode();
    void reshuffleLayout(bool hidePlots);
};

#endif // MAINWINDOW_H
