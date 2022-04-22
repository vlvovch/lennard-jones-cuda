#include <QtGui>
#include "mainwindow.h"
#include "cuda_runtime.h"

#include <algorithm>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    QWidget *centW = new QWidget;
    setCentralWidget(centW);

    QGLFormat fmt;
    fmt.setAccum(true);

    glWidget = new GLWidget(fmt);
    connect(glWidget, SIGNAL(updateGraph()), this, SLOT(updateEner()));

    //QVBoxLayout *
    lhsLayout = new QVBoxLayout;
    //QGridLayout *
    plotsLayout = new QGridLayout;

    ploten = new QCustomPlot;
    ploten->xAxis->setLabel("t");
    ploten->yAxis->setLabel("");
    ploten->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    plotrdf = new QCustomPlot;
    plotrdf->xAxis->setLabel("r/sigma");
    plotrdf->yAxis->setLabel(tr("g(r)"));
    plotrdf->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    plotvel = new QCustomPlot;
    plotvel->xAxis->setLabel("v");
    plotvel->yAxis->setLabel("f(v)");
    plotvel->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    plotsLayout->addWidget(glWidget, 0, 0);
    plotsLayout->addWidget(plotrdf, 0, 1);
    plotsLayout->addWidget(ploten, 1, 0);
    plotsLayout->addWidget(plotvel, 1, 1);
    //glWidget->setFormat(fmt);

    QHBoxLayout *layCheckPlots = new QHBoxLayout;
    QCheckBox *CBhideplots = new QCheckBox(tr("Maximize 3D view"));
    CBhideplots->setChecked(false);
    connect(CBhideplots, SIGNAL(toggled(bool)), this, SLOT(reshuffleLayout(bool)));
    layCheckPlots->addWidget(CBhideplots, 0, Qt::AlignRight);

    plotsWidget = new QWidget;
    plotsWidget->setLayout(plotsLayout);
    //lhsLayout->addLayout(plotsLayout);
    lhsLayout->addWidget(plotsWidget);
    lhsLayout->addLayout(layCheckPlots);


    QGroupBox *ParamGroup = new QGroupBox(tr("System parameters"));
    QGridLayout *paramsLayout = new QGridLayout;
    paramsLayout->setAlignment(Qt::AlignLeft);

    QLabel *labelN = new QLabel(tr("N = "));
    spinN = new QSpinBox();
    spinN->setRange(1, 100000);
    spinN->setValue(100);
    paramsLayout->addWidget(labelN,0,0);
    paramsLayout->addWidget(spinN,0,1);

    QLabel* labelT = new QLabel(tr("T<sup>*</sup> = "));
    spinT = new QDoubleSpinBox();
    spinT->setDecimals(3);
    spinT->setRange(0., 1000.);
    //spinT->setValue(0.722);
    spinT->setValue(0.710);
    spinT->setSingleStep(0.1);
    paramsLayout->addWidget(labelT, 1, 0);
    paramsLayout->addWidget(spinT, 1, 1);

    QLabel *labelrho = new QLabel(tr("ρ<sup>*</sup> = "));
    spinrho = new QDoubleSpinBox();
    spinrho->setDecimals(3);
    spinrho->setRange(0., 1000.);
    spinrho->setValue(0.844);
    spinrho->setSingleStep(0.1);
    paramsLayout->addWidget(labelrho,2,0);
    paramsLayout->addWidget(spinrho,2,1);

    ParamGroup->setLayout(paramsLayout);

    QGroupBox *RenderGroup = new QGroupBox(tr("Render method"));
    QHBoxLayout *renderLayout = new QHBoxLayout;
    RBSpheres = new QRadioButton(tr("Spheres"));
    RBPoints = new QRadioButton(tr("Points"));
    renderLayout->addWidget(RBSpheres);
    renderLayout->addWidget(RBPoints);
    RenderGroup->setLayout(renderLayout);
    RBSpheres->setChecked(true);
    connect(RBSpheres, SIGNAL(toggled(bool)), this, SLOT(updateRenderMode()));

    QGroupBox *EnsGroup = new QGroupBox(tr("Ensemble"));
    QHBoxLayout *ensLayout = new QHBoxLayout;
    RBMicro = new QRadioButton(tr("Microcanonical"));
    RBCanon = new QRadioButton(tr("Canonical"));
    ensLayout->addWidget(RBMicro);
    ensLayout->addWidget(RBCanon);
    EnsGroup->setLayout(ensLayout);
    RBMicro->setChecked(true);
    connect(RBCanon, SIGNAL(toggled(bool)), this, SLOT(updateEnsemble()));

    QGroupBox* BoundaryGroup = new QGroupBox(tr("Boundary conditions"));
    QHBoxLayout* bndLayout = new QHBoxLayout;
    RBPeriodic = new QRadioButton(tr("Periodic"));
    RBPeriodic->setToolTip(tr("Minimum image convention"));
    RBWalls = new QRadioButton(tr("Hard walls"));
    RBExpansion = new QRadioButton(tr("Expansion"));
    bndLayout->addWidget(RBPeriodic);
    bndLayout->addWidget(RBWalls);
    bndLayout->addWidget(RBExpansion);
    BoundaryGroup->setLayout(bndLayout);
    RBPeriodic->setChecked(true);
    connect(RBPeriodic, SIGNAL(toggled(bool)), this, SLOT(updateBoundary()));
    connect(RBWalls, SIGNAL(toggled(bool)), this, SLOT(updateBoundary()));
    connect(RBExpansion, SIGNAL(toggled(bool)), this, SLOT(updateBoundary()));

    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
	  printf("Number of CUDA-supported devices: %d\n", deviceCount);

    QGroupBox *HardGroup = new QGroupBox(tr("Calculate on"));
    QHBoxLayout *hardLayout = new QHBoxLayout;
    RBCPU = new QRadioButton(tr("CPU"));
    RBCUDA = new QRadioButton(tr("GPU (CUDA)"));
    if (deviceCount==0) RBCUDA->setEnabled(false);
    hardLayout->addWidget(RBCPU);
    hardLayout->addWidget(RBCUDA);
    HardGroup->setLayout(hardLayout);
    RBCPU->setChecked(true);
    connect(RBCUDA, SIGNAL(toggled(bool)), this, SLOT(updateHardwareMode()));

    QGroupBox *SimulGroup = new QGroupBox(tr("Simulation parameters"));
    QGridLayout *simulGridLayout = new QGridLayout;
    simulGridLayout->setAlignment(Qt::AlignLeft);

    QLabel *labeldtau = new QLabel(tr("dt<sup>*</sup> = "));
    spindtau = new QDoubleSpinBox();
    spindtau->setDecimals(4);
    spindtau->setRange(0., 1000.);
    spindtau->setValue(0.004);
    spindtau->setSingleStep(0.001);
    simulGridLayout->addWidget(labeldtau,0,0);
    simulGridLayout->addWidget(spindtau,0,1);

    QLabel *labelpre = new QLabel(tr("Iters per frame = "));
    spinpre = new QSpinBox();
    spinpre->setRange(1, 100000);
    spinpre->setValue(10);
    simulGridLayout->addWidget(labelpre,1,0);
    simulGridLayout->addWidget(spinpre,1,1);
    SimulGroup->setLayout(simulGridLayout);


    QVBoxLayout *buttonLayout = new QVBoxLayout;
    buttonLayout->setAlignment(Qt::AlignLeft);
    buttonCont = new QPushButton(tr("Continue"));
    buttonStop = new QPushButton(tr("Pause"));
    buttonRestart = new QPushButton(tr("Restart"));
    buttonReaverage = new QPushButton(tr("Restart Averaging"));
    //buttonRestart->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    buttonStop->setEnabled(false);
    connect(buttonCont, SIGNAL(clicked()), this, SLOT(cont()));
    connect(buttonStop, SIGNAL(clicked()), this, SLOT(stop()));
    connect(buttonRestart, SIGNAL(clicked()), this, SLOT(restart()));
    connect(buttonReaverage, SIGNAL(clicked()), this, SLOT(reav()));
    buttonLayout->addWidget(buttonCont);
    buttonLayout->addWidget(buttonStop);
    buttonLayout->addWidget(buttonRestart);
    buttonLayout->addWidget(buttonReaverage);
    QLabel *labelCopyright = new QLabel(tr("© 2013-2021 V. Vovchenko"));


    QVBoxLayout *interLayout = new QVBoxLayout;
    //interLayout->setAlignment(Qt::AlignTop);
    interLayout->addWidget(ParamGroup);
    interLayout->addWidget(RenderGroup);
    interLayout->addWidget(EnsGroup);
    interLayout->addWidget(BoundaryGroup);
    //interLayout->addWidget(IntGroup);
    interLayout->addWidget(HardGroup);
    interLayout->addWidget(SimulGroup);
    //interLayout->addWidget(buttonRestart);
    interLayout->addLayout(buttonLayout);
    //interLayout->addLayout(traceGridLayout);
    interLayout->addWidget(labelCopyright, 0, Qt::AlignRight);

    QWidget *interWidget = new QWidget();
    interWidget->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    interWidget->setLayout(interLayout);


    QHBoxLayout *mainLayout = new QHBoxLayout;
    mainLayout->addLayout(lhsLayout);
    mainLayout->addWidget(interWidget, 0, Qt::AlignTop);
    centralWidget()->setLayout(mainLayout);

    glWidget->dt = spindtau->value();
    glWidget->pre = spinpre->value();
    glWidget->systGL = new MDSystemGL(getConfig(), RBSpheres->isChecked());

    ploten->xAxis->setRange(0., 200.);
    ploten->yAxis->setRange(0., glWidget->systGL->m_system->m_config.T0 * 2);
    ploten->addGraph();
    ploten->graph(0)->setName(tr("Energy per particle"));
    ploten->graph(0)->setPen(QPen(Qt::black, 4));
    ploten->graph(0)->addData(0., (glWidget->systGL->m_system->K+glWidget->systGL->m_system->V)/glWidget->systGL->m_system->m_config.N);
    ploten->addGraph();
    ploten->graph(1)->setName(tr("Instant temperature"));
    ploten->graph(1)->setPen(QPen(Qt::red, 4));
    ploten->graph(1)->addData(0., glWidget->systGL->m_system->T);
    ploten->legend->setVisible(true);
    //ploten->addGraph();
    //ploten->graph(2)->setName(tr("Compressibility"));
    //ploten->graph(2)->setPen(QPen(Qt::blue, 4));
    //ploten->graph(2)->addData(0., glWidget->systGL->m_system->P / (glWidget->systGL->m_system->m_config.T0 * glWidget->systGL->m_system->m_config.rho));
    //ploten->legend->setVisible(true);

    plotrdf->yAxis->setRange(0., 3.2);
    plotrdf->addGraph();
    plotrdf->graph(0)->setName(tr("Radial distribution function"));
    plotrdf->graph(0)->setPen(QPen(Qt::black, 4));
    plotrdf->addGraph();
    plotrdf->graph(1)->setName(tr("Low-density limit, g(r) = exp[-u(r)/T]"));
    plotrdf->graph(1)->setPen(QPen(Qt::red, 4, Qt::DashLine));
    plotrdf->legend->setVisible(true);
    rdfiter = 0;
    //drawRDF();
    QCPItemLine* item = new QCPItemLine(plotrdf);
    plotrdf->addItem(item);
    item->setPen(QPen(Qt::black, 4, Qt::DashLine));
    item->start->setCoords(0, 1);
    item->end->setCoords(5, 1);


    plotvel->yAxis->setRange(0., 0.8);
    plotvel->xAxis->setRange(0., 12.);
    plotvel->addGraph();
    plotvel->graph(0)->setName(tr("Velocity distribution"));
    plotvel->graph(0)->setPen(QPen(Qt::black, 4));
    plotvel->addGraph();
    plotvel->graph(1)->setName(tr("Maxwell-Boltzmann"));
    plotvel->graph(1)->setPen(QPen(Qt::red, 4));
    plotvel->legend->setVisible(true);
    drawvel();
    drawRDF();

    setWindowTitle(tr("Lennard-Jones Molecular Dynamics"));
}

MainWindow::~MainWindow()
{

}

MDSystem::MDSystemConfiguration MainWindow::getConfig()
{
    MDSystem::MDSystemConfiguration config;

    config.N = spinN->value();
    config.T0 = spinT->value();
    config.rho = spinrho->value();
    config.canonical = RBCanon->isChecked();
    config.useCUDA = RBCUDA->isChecked();
    config.CUDABlockSize = 256;
    //config.boundaryConditions = !RBPeriodic->isChecked();
    config.boundaryConditions = getBoundaryCondition();

    return config;
}

void MainWindow::keyPressEvent(QKeyEvent *e)
{
    if (e->key() == Qt::Key_Escape)
        close();
    else
        QWidget::keyPressEvent(e);
}

void MainWindow::resizeEvent(QResizeEvent *event)
{

}

void MainWindow::updateEner()
{
    ploten->graph(0)->addData(glWidget->systGL->m_system->getTime(), (glWidget->systGL->m_system->K+glWidget->systGL->m_system->V)/glWidget->systGL->m_system->m_config.N);
    ploten->graph(1)->addData(glWidget->systGL->m_system->getTime(), glWidget->systGL->m_system->T);
    //ploten->graph(2)->addData(glWidget->systGL->m_system->getTime(), glWidget->systGL->m_system->P / (glWidget->systGL->m_system->m_config.T0 * glWidget->systGL->m_system->m_config.rho));
    if (glWidget->systGL->m_system->getTime()>200.) 
      ploten->xAxis->setRange(ploten->xAxis->range().lower+glWidget->dt, ploten->xAxis->range().upper+glWidget->dt);
    

    double ymin = ploten->yAxis->range().lower;
    double ymax = ploten->yAxis->range().upper;

    double T = glWidget->systGL->m_system->T;
    double u = (glWidget->systGL->m_system->K + glWidget->systGL->m_system->V) / glWidget->systGL->m_system->m_config.N;
    double Z = 1.;// glWidget->systGL->m_system->P / (glWidget->systGL->m_system->m_config.T0 * glWidget->systGL->m_system->m_config.rho);
    ymin = std::min(ymin, std::min(T, Z));
    if (u < 0 && ymin > u)
      ymin = u * 1.2;
    ymax = std::max(ymax, std::max(T, Z));
    ymax = std::max(ymax, u);


    //ploten->yAxis->setRange(0., glWidget->systGL->m_system->m_config.T0 * 2);
    ploten->yAxis->setRange(ymin, ymax);

    ploten->replot();

    drawRDF();
    plotrdf->replot();
    drawvel();
    plotvel->replot();
}

void MainWindow::updateRenderMode()
{
  glWidget->systGL->setRenderMode(RBSpheres->isChecked());
}

void MainWindow::updateEnsemble()
{
  glWidget->systGL->setCanonical(RBCanon->isChecked());
}

void MainWindow::updateBoundary()
{
  glWidget->systGL->setBoundaryCondition(getBoundaryCondition());
}

void MainWindow::updateHardwareMode()
{
  glWidget->systGL->setHardwareMode(RBCUDA->isChecked());
}

void MainWindow::drawRDF()
{
    plotrdf->graph(0)->clearData();
    SplineFunction rdf = glWidget->systGL->m_system->RDF(plotrdf->xAxis->range().upper, 0.03);
    if (rdfiter==0)
    {
        rdftot = rdf;
    }
    else
    {
        for(int i=0;i<rdf.vals.size();++i)
            rdftot.vals[i].second = (rdftot.vals[i].second*rdfiter + rdf.vals[i].second) / (rdfiter+1);
    }
    rdfiter++;
    plotrdf->graph(0)->addData(0., 0.);
    double maxy = 0.;
    for (int i = 0; i < rdftot.vals.size(); ++i) {
      plotrdf->graph(0)->addData(rdftot.vals[i].first, rdftot.vals[i].second);
      if (rdftot.vals.size() - 1 != i)
        maxy = std::max(rdftot.vals[i].second, maxy);
    }

    if (maxy > 3.2) {
      plotrdf->yAxis->setRange(0., maxy * 1.2);
    }
    else {
      plotrdf->yAxis->setRange(0., 3.2);
    }

    plotrdf->graph(1)->clearData();
    int iters = 100;
    double dr = plotrdf->xAxis->range().upper / iters;
    for (int i = 0; i < iters; ++i)
    {
      double r = dr * (i + 0.5);
      double U = 4. * (1. / pow(r, 12) - 1. / pow(r, 6));
      double T = glWidget->systGL->m_system->T;
      double g = exp(-U / T);
      plotrdf->graph(1)->addData(r, g);
    }
}

void MainWindow::drawvel()
{
    plotvel->graph(0)->clearData();
    plotvel->graph(1)->clearData();
    glWidget->systGL->m_system->updatevelo();
    SplineFunction velf = glWidget->systGL->m_system->getvelo();
    for(int i=0;i<velf.vals.size();++i)
        plotvel->graph(0)->addData(velf.vals[i].first, velf.vals[i].second);
    int iters = 101;
    double dv = plotvel->xAxis->range().upper / (iters-1);
    for(int i=0;i<iters;++i)
    {
        plotvel->graph(1)->addData(dv*i, glWidget->systGL->m_system->Maxwell(dv*i));
    }
}


void MainWindow::cont()
{
    buttonCont->setEnabled(false);
    buttonStop->setEnabled(true);
    glWidget->anim->start(10);
    RBPeriodic->setEnabled(false);
    RBWalls->setEnabled(false);
    RBExpansion->setEnabled(true);
}

void MainWindow::stop()
{
    buttonCont->setEnabled(true);
    buttonStop->setEnabled(false);
    glWidget->anim->stop();
}

void MainWindow::reav()
{
    glWidget->systGL->m_system->initvelo();
    glWidget->systGL->m_system->resetAveraging();
    rdfiter = 0;
}

void MainWindow::restart()
{
    glWidget->dt = spindtau->value();
    glWidget->pre = spinpre->value();
    glWidget->systGL->Reinitialize(getConfig(), RBSpheres->isChecked());
    ploten->xAxis->setRange(0., 200.);
    ploten->yAxis->setRange(0., glWidget->systGL->m_system->m_config.T0 * 2.);
    ploten->graph(0)->clearData();
    ploten->graph(0)->addData(0., (glWidget->systGL->m_system->K+glWidget->systGL->m_system->V)/glWidget->systGL->m_system->m_config.N);
    ploten->graph(1)->clearData();
    ploten->graph(1)->addData(0., glWidget->systGL->m_system->T);
    ploten->replot();
    rdfiter = 0;
    plotrdf->replot();
    drawvel();
    plotvel->replot();
    glWidget->updateGL();
    RBPeriodic->setEnabled(true);
    RBWalls->setEnabled(true);
}

void MainWindow::reshuffleLayout(bool hidePlots) {
  if (hidePlots) {
    lhsLayout->removeWidget(plotsWidget);
    lhsLayout->insertWidget(0, glWidget);
  }
  else {
    lhsLayout->removeWidget(glWidget);
    plotsLayout->addWidget(glWidget, 0, 0);
    lhsLayout->insertWidget(0, plotsWidget);
  }
}

int MainWindow::getBoundaryCondition()
{
  if (RBPeriodic->isChecked())
    return 0;
  if (RBWalls->isChecked())
    return 1;
  return 2;
}
