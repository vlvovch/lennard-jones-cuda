#include "MDSystem.h"
#include <cmath>
#include <ctime>
#include <cstring>
#include <fstream>
#include <iostream>
#define PI 3.141592653589793238462643

extern "C"
{
    //void checkCUDA();
    void allocateNBodyArrays(float* vel[2], int numBodies);
    void allocateArray(float** dest, int number);
    void deleteNBodyArrays(float* vel[2]);
    void deleteArray(float* arr);
    void calculateNForces(float* Pos, float* Force, float *host_pressure,
      int numBodies, float host_L, int Lperiodic, int* host_RDF, float host_dr2, int p, int q);
    void copyArrayFromDevice(float* host, const float* device, unsigned int pbo, int numBodies);
    void copyArrayToDevice(float* device, const float* host, int numBodies);
    void registerGLBufferObject(unsigned int pbo);
    void unregisterGLBufferObject(unsigned int pbo);
    void threadSync();
    void threadExit();
}


double MDSystem::MaxwellDistributionGenerator::prob(double alpha)
{
    double tlog = -log(alpha);
    return tlog * tlog * exp(-tlog*tlog) / alpha;
}

double MDSystem::MaxwellDistributionGenerator::sampleAlphaV()
{
    bool fl = false;
    double talpha = 0.;
    while (!fl) {
        talpha = rangen.randDblExc();
        double eta2 = rangen.rand() * probmax;
        if (eta2 < prob(talpha))
            fl = true;
    }
    return talpha;
}

MDSystem::MDSystem(const MDSystem::MDSystemConfiguration& config)
{
    CUDAInit = false;

    h_Pos = NULL;
    h_Vel = NULL;
    h_Force = NULL;

    Reinitialize(config);
}

void MDSystem::Reinitialize(const MDSystem::MDSystemConfiguration &config)
{
    m_config = config;

    L = pow(m_config.N / m_config.rho, 1./3.);

    this->T = m_config.T0;

    t = 0;


    if (h_Pos != NULL)
        delete [] h_Pos;
    if (h_Vel != NULL)
        delete [] h_Vel;
    if (h_Force != NULL)
        delete [] h_Force;

    const int &N = m_config.N;

    h_Pos = new float[4*N];
    h_Vel = new float[4*N];
    h_Force = new float[4*N];

    SampleInitialConditions();

    // RDF
    rdf_dr2 = std::max(0.2 * sqrt(100. / N), 0.05);
    if (250 * rdf_dr2 < 25.0)
      rdf_dr2 = 25.0 / 250;

    NdNdr2 = std::vector<int>(256, 0);

    ReallocateMemory();

    CalculateForces();

    CalculateParameters();

    veloIters = 0;

    initvelo();

    momN = momN2 = momN3 = momN4 = 0.;
    momiters = 0;

    av_U_tot = av_T_tot = av_p_tot = 0.;
    av_iters = 0;
}

void MDSystem::ReallocateMemory()
{
  const int& N = m_config.N;
  if (m_config.useCUDA)
  {
    if (CUDAInit)
    {
      deleteArray(d_Pos);
      deleteArray(d_Vel);
      deleteArray(d_Force);
    }
    allocateArray(&d_Pos, N);
    copyArrayToDevice(d_Pos, h_Pos, N);
    allocateArray(&d_Vel, N);
    copyArrayToDevice(d_Vel, h_Vel, N);
    allocateArray(&d_Force, N);
    CUDAInit = true;
  }
  else
  {
    if (CUDAInit)
    {
      deleteArray(d_Pos);
      deleteArray(d_Vel);
      deleteArray(d_Force);
    }
    CUDAInit = false;
  }
}


void MDSystem::SampleInitialConditions()
{
  const int& N = m_config.N;
  int Nsingle = ceil(pow(N, 1. / 3.));
  double dL = L / Nsingle;
  bool fl = 0;
  std::vector<int> indis;
  for (int i = 0; i < 4 * N; i += 4)
    indis.push_back(i);
  std::random_shuffle(indis.begin(), indis.end());
  for (int ii = 0; ii < N; ii++)
  {
    int i = indis[ii];
    int iN = i / 4;
    int ix = iN % Nsingle;
    int iy = (iN / Nsingle) % Nsingle;
    int iz = iN / (Nsingle * Nsingle);

    h_Pos[i] = (ix + 0.5) * dL;
    h_Pos[i + 1] = (iy + 0.5) * dL;
    h_Pos[i + 2] = (iz + 0.5) * dL;
    h_Pos[i + 3] = L / 150.f;
    

    double tmpv, tmpvcth, tmpvph;
    tmpv = sqrt(2. * T) * (-log(m_MaxwellGenerator.sampleAlphaV()));
    tmpvcth = 2. * rangen.rand() - 1.;
    tmpvph = 2 * PI * rangen.rand();
    h_Vel[i] = tmpv * sqrt(1 - tmpvcth * tmpvcth) * cos(tmpvph);
    h_Vel[i + 1] = tmpv * sqrt(1 - tmpvcth * tmpvcth) * sin(tmpvph);
    h_Vel[i + 2] = tmpv * tmpvcth;
  }
  CorrectTotalMomentum();
  RenormalizeVelocities(true);
}

void MDSystem::CorrectTotalMomentum()
{
  double tPx = 0., tPy = 0., tPz = 0., totmass = 0.;
  for (int i = 0; i < 4 * m_config.N; i += 4)
  {
    tPx += h_Vel[i];
    tPy += h_Vel[i + 1];
    tPz += h_Vel[i + 2];
    totmass += 1.;
  }

  //for debugging
  //printf("%lf %lf %lf\n", tPx, tPy, tPz);

  for (int i = 0; i < 4 * m_config.N; i += 4)
  {
    h_Vel[i] += -tPx / totmass;
    h_Vel[i + 1] += -tPy / totmass;
    h_Vel[i + 2] += -tPz / totmass;
  }

  tPx = tPy = tPz = 0.;

  for (int i = 0; i < 4 * m_config.N; i += 4)
  {
    tPx += h_Vel[i];
    tPy += h_Vel[i + 1];
    tPz += h_Vel[i + 2];
    totmass += 1.;
  }

  //for debugging
  //printf("%lf %lf %lf\n", tPx, tPy, tPz);
}

MDSystem::~MDSystem(void)
{
    delete [] h_Pos;
    delete [] h_Vel;
    delete [] h_Force;
    if (CUDAInit)
    {
        deleteArray(d_Pos);
        deleteArray(d_Vel);
        deleteArray(d_Force);
    }
    threadExit();
}

void MDSystem::CalculateForces()
{
  CalculateForces(h_Force);
}

void MDSystem::CalculateForces(float *frc)
{
    const int &N = m_config.N;
    if (m_config.useCUDA)
    {
        copyArrayToDevice(d_Pos, h_Pos, N);

        float pressure = 0.f;
        calculateNForces(d_Pos, d_Force, &pressure,
          N, static_cast<float>(L), !m_config.boundaryConditions, &NdNdr2[0], rdf_dr2, m_config.CUDABlockSize, 1);

        P = pressure;

        copyArrayFromDevice(frc, d_Force, 0, N);
    }
    else
    {
        std::fill(NdNdr2.begin(), NdNdr2.end(), 0);
      
        double r2,r6;
        V = 0.;
        P = 0.;
        Pshear = 0.;
        for(int i=0;i<4*N;i+=4)
        {
            frc[i] = 0.;
            frc[i+1] = 0.;
            frc[i+2] = 0.;
            for(int j=0;j<4*N;j+=4)
            {
                if (j!=i)
                {
                  float rx = h_Pos[i] - h_Pos[j];
                  float ry = h_Pos[i + 1] - h_Pos[j + 1];
                  float rz = h_Pos[i + 2] - h_Pos[j + 2];

                  if (m_config.boundaryConditions == 0) {
                    rx = rx - L * fast_round(rx / L);
					          ry = ry - L * fast_round(ry / L);
					          rz = rz - L * fast_round(rz / L);
                  }

                  r2 = rx * rx + ry * ry + rz * rz;

                  {
                    int indrdf = floor(r2 / rdf_dr2);
                    if (indrdf < NdNdr2.size())
                      NdNdr2[indrdf]++;
                  }

                  //if (r2<rc*rc)
                  {
                      r2 = 1/r2;
                      r6 = r2*r2*r2;
                      double fijx = r2 * (12 * r6 * r6 - 6 * r6) * rx;
                      double fijy = r2 * (12 * r6 * r6 - 6 * r6) * ry;
                      double fijz = r2 * (12 * r6 * r6 - 6 * r6) * rz;
                      frc[i] += fijx;
                      frc[i + 1] += fijy;
                      frc[i + 2] += fijz;
                      V += (1*r6*r6 - 1*r6);// - Vcut;
                      P += rx * fijx + ry * fijy + rz * fijz;
                      Pshear += -rx * fijy;
                  }
                }
            }
            frc[i] *= 4;
            frc[i+1] *= 4;
            frc[i+2] *= 4;
        }
        P *= 4. / 3. / 2.;
        V *= 4. / 2.;
        Pshear *= 4. / 2.;
    }
}

double MDSystem::CalculateXi(float* frc, float* vel)
{
  const int& N = m_config.N;
  double retnum = 0.f, retzn = 0.f;
  for (int i = 0; i < 4 * N; i += 4)
  {
    retnum += vel[i] * frc[i] + vel[i + 1] * frc[i + 1] + vel[i + 2] * frc[i + 2];
    retzn += vel[i] * vel[i] + vel[i + 1] * vel[i + 1] + vel[i + 2] * vel[i + 2];
  }
  return retnum / retzn;
}

void MDSystem::CalculateParameters()
{
    const int &N = m_config.N;
    K = 0.;
    
    if (!m_config.useCUDA)
    {
        for(int i=0;i<4*N;i+=4)
        {
            K += (h_Vel[i]*h_Vel[i] + h_Vel[i+1]*h_Vel[i+1] + h_Vel[i+2]*h_Vel[i+2]) / 2.;
            Pshear += -h_Vel[i] * h_Vel[i+1];
        }
    }
    else
    {
        V = 0.;
        for(int i=0;i<4*N;i+=4)
        {
            K += (h_Vel[i]*h_Vel[i] + h_Vel[i+1]*h_Vel[i+1] + h_Vel[i+2]*h_Vel[i+2]) / 2.;
            V += h_Force[i+3];
        }
        V *= 4 / 2;
    }
    T = 2. * K / 3. / N;
    P += m_config.N * T;
    P /= (m_config.N / m_config.rho);
    U = K + V;

    Pshear /= (m_config.N / m_config.rho);

    av_iters++;
    av_U_tot += U;
    av_p_tot += P;
    av_T_tot += T;
}

double MDSystem::KineticTemperature(float* vel)
{
  double ret = 0.;

  for (int i = 0; i < 4 * m_config.N; i += 4)
  {
    ret += (vel[i] * vel[i] + vel[i + 1] * vel[i + 1] + vel[i + 2] * vel[i + 2]);
  }

  ret *= 1. / 3. / m_config.N;

  return ret;
}

void MDSystem::RenormalizeVelocities(bool RecalculateTkin)
{
  double Tkin = T;
  if (RecalculateTkin)
    Tkin = KineticTemperature(h_Vel);
  for(int i=0;i<4*m_config.N;i+=4)
  {
      h_Vel[i] *= sqrt(m_config.T0/ Tkin);
      h_Vel[i+1] *= sqrt(m_config.T0/ Tkin);
      h_Vel[i+2] *= sqrt(m_config.T0/ Tkin);
  }
  K *= m_config.T0/ Tkin;
  T = m_config.T0;
  U = K + V;
}

void MDSystem::RenormalizeVelocitiesToEnergy(double ust)
{
  double Kold = K;
  double Udes = ust * m_config.N;
  double Kdes = Udes - V;
  for (int i = 0; i < 4 * m_config.N; i += 4)
  {
    h_Vel[i] *= sqrt(Kdes / Kold);
    h_Vel[i + 1] *= sqrt(Kdes / Kold);
    h_Vel[i + 2] *= sqrt(Kdes / Kold);
  }
  K = Kdes;
  U = K + V;
}

void MDSystem::ApplyBoundaryConditions()
{
  
  if (m_config.boundaryConditions == 2) {
    // No boundary conditions (expansion)
    return;
  }
  
  if (m_config.boundaryConditions == 0) {
    for (int i = 0; i < 4 * m_config.N; i += 4)
    {
      if (h_Pos[i]     < 0.)  h_Pos[i] += L;
      if (h_Pos[i]     > L )  h_Pos[i] -= L;
      if (h_Pos[i + 1] < 0.)  h_Pos[i + 1] += L;
      if (h_Pos[i + 1] > L )  h_Pos[i + 1] -= L;
      if (h_Pos[i + 2] < 0.)  h_Pos[i + 2] += L;
      if (h_Pos[i + 2] > L )  h_Pos[i + 2] -= L;
    }
  }
  else {
    for (int i = 0; i < 4 * m_config.N; i += 4)
    {
      if (h_Pos[i]     < 0. && h_Vel[i] < 0)   h_Vel[i] = -h_Vel[i];
      if (h_Pos[i]     > L  && h_Vel[i] > 0)   h_Vel[i] = -h_Vel[i];
      if (h_Pos[i + 1] < 0. && h_Vel[i + 1] < 0) h_Vel[i + 1] = -h_Vel[i + 1];
      if (h_Pos[i + 1] > L  && h_Vel[i + 1] > 0) h_Vel[i + 1] = -h_Vel[i + 1];
      if (h_Pos[i + 2] < 0. && h_Vel[i + 2] < 0) h_Vel[i + 2] = -h_Vel[i + 2];
      if (h_Pos[i + 2] > L  && h_Vel[i + 2] > 0) h_Vel[i + 2] = -h_Vel[i + 2];
    }
  }
}

void MDSystem::Integrate(double dt)
{
    const int &N = m_config.N;

    if (!m_config.canonical) {
      for (int i = 0; i < 4 * N; i += 4)
      {
        float fx = h_Force[i], fy = h_Force[i + 1], fz = h_Force[i + 2];

        h_Pos[i] += dt * h_Vel[i] + dt * dt * fx /*/ m *// 2.;
        h_Pos[i + 1] += dt * h_Vel[i + 1] + dt * dt * fy /*/ m*/ / 2.;
        h_Pos[i + 2] += dt * h_Vel[i + 2] + dt * dt * fz /*/ m*/ / 2.;

        h_Vel[i] += dt * fx /*/ m*/ / 2.;
        h_Vel[i + 1] += dt * fy /*/ m*/ / 2.;
        h_Vel[i + 2] += dt * fz /*/ m*/ / 2.;
      }

      CalculateForces();

      for (int i = 0; i < 4 * N; i += 4)
      {
        h_Vel[i] += dt * h_Force[i] /*/ m*/ / 2.;
        h_Vel[i + 1] += dt * h_Force[i + 1] /*/ m*/ / 2.;
        h_Vel[i + 2] += dt * h_Force[i + 2] /*/ m*/ / 2.;
      }
    }
    else if (1) {
      float* t_Force = new float[4 * N];
      float* t_Vel = new float[4 * N];

      for (int i = 0; i < 4 * N; i += 4)
      {
        float fx = h_Force[i], fy = h_Force[i + 1], fz = h_Force[i + 2];

        h_Pos[i] += dt * h_Vel[i] + dt * dt * fx /*/ m*/ / 2.;
        h_Pos[i + 1] += dt * h_Vel[i + 1] + dt * dt * fy /*/ m*/ / 2.;
        h_Pos[i + 2] += dt * h_Vel[i + 2] + dt * dt * fz /*/ m*/ / 2.;

        t_Force[i] = 0.5f * h_Force[i];
        t_Force[i+1] = 0.5f * h_Force[i + 1];
        t_Force[i+2] = 0.5f * h_Force[i + 2];
      }

      CalculateForces();

      for (int i = 0; i < 4 * N; i += 4)
      {
        t_Force[i] += 0.5f * h_Force[i];
        t_Force[i + 1] += 0.5f * h_Force[i + 1];
        t_Force[i + 2] += 0.5f * h_Force[i + 2];
      }

      for (int i = 0; i < 4 * N; i += 4)
      {
        t_Vel[i] = h_Vel[i] + dt * t_Force[i] /*/ m*/ / 2.;
        t_Vel[i + 1] = h_Vel[i + 1] + dt * t_Force[i+1] /*/ m*/ / 2.;
        t_Vel[i + 2] = h_Vel[i + 2] + dt * t_Force[i+2] /*/ m*/ / 2.;
      }

      double Tkin = KineticTemperature(t_Vel);
      double chi = sqrt(m_config.T0 / Tkin);

      for (int i = 0; i < 4 * N; i += 4)
      {
        h_Vel[i] = (2. * chi - 1.) * h_Vel[i] + chi * dt * t_Force[i] /*/ m*/ ;
        h_Vel[i + 1] = (2. * chi - 1.) * h_Vel[i + 1] + chi * dt * t_Force[i + 1] /*/ m*/;
        h_Vel[i + 2] = (2. * chi - 1.) * h_Vel[i + 2] + chi * dt * t_Force[i + 2] /*/ m*/;
      }

      delete [] t_Force;
      delete [] t_Vel;
    }
    else { // alternative canonical ensemble, not used anymore
      float* t_Force = new float[4 * N];
      float* t_Vel = new float[4 * N];

      double txi = CalculateXi(h_Force, h_Vel);

      for (int i = 0; i < 4 * N; i += 4)
      {
        t_Vel[i] = h_Vel[i] - h_Vel[i] * txi * dt + dt * h_Force[i] /*/ m*/;
        t_Vel[i+1] = h_Vel[i+1] - h_Vel[i+1] * txi * dt + dt * h_Force[i+1] /*/ m*/;
        t_Vel[i+2] = h_Vel[i+2] - h_Vel[i+2] * txi * dt + dt * h_Force[i+2] /*/ m*/;
        
        h_Vel[i] += -h_Vel[i] * txi * dt / 2. + dt * h_Force[i] / 2. /*/ m*/;
        h_Vel[i+1] += -h_Vel[i+1] * txi * dt / 2. + dt * h_Force[i+1] / 2. /*/ m*/;
        h_Vel[i+2] += -h_Vel[i+2] * txi * dt / 2. + dt * h_Force[i+2] / 2. /*/ m*/;

        h_Pos[i] += dt * h_Vel[i];
        h_Pos[i + 1] += dt * h_Vel[i + 1];
        h_Pos[i + 2] += dt * h_Vel[i + 2];

        
        
        //float fx = h_Force[i], fy = h_Force[i + 1], fz = h_Force[i + 2];

        //fx += -txi * m * h_Vel[i];
        //fy += -txi * m * h_Vel[i + 1];
        //fz += -txi * m * h_Vel[i + 2];

        //h_Pos[i] += dt * h_Vel[i] + dt * dt * fx / m / 2.;
        //h_Pos[i + 1] += dt * h_Vel[i + 1] + dt * dt * fy / m / 2.;
        //h_Pos[i + 2] += dt * h_Vel[i + 2] + dt * dt * fz / m / 2.;

        //h_Vel[i] += dt * fx / m;
        //h_Vel[i + 1] += dt * fy / m;
        //h_Vel[i + 2] += dt * fz / m;

        //t_Force[i] = fx;
        //t_Force[i + 1] = fy;
        //t_Force[i + 2] = fz;
      }

      CalculateForces();

      //double Tkin = KineticTemperature(t_Vel);
      txi = CalculateXi(h_Force, t_Vel);
      //txi = (sqrt(Tkin / m_config.T0) - 1.) / (dt / 2.);

      for (int i = 0; i < 4 * N; i += 4)
      {
        h_Vel[i] = (h_Vel[i] + dt / 2. * h_Force[i] /*/ m*/) /  (1. + txi * dt / 2.);
        h_Vel[i + 1] = (h_Vel[i + 1] + dt / 2. * h_Force[i + 1] /*/ m*/) / (1. + txi * dt / 2.);
        h_Vel[i + 2] = (h_Vel[i + 2] + dt / 2. * h_Force[i + 2] /*/ m*/) / (1. + txi * dt / 2.);
        //float fx = h_Force[i], fy = h_Force[i + 1], fz = h_Force[i + 2];

        //fx += -txi * m * h_Vel[i];
        //fy += -txi * m * h_Vel[i + 1];
        //fz += -txi * m * h_Vel[i + 2];

        //h_Vel[i] += dt * (fx - t_Force[i]) / m / 2.;
        //h_Vel[i + 1] += dt * (fy - t_Force[i + 1]) / m / 2.;
        //h_Vel[i + 2] += dt * (fz - t_Force[i + 2]) / m / 2.;
      }

      delete[] t_Force;
      delete[] t_Vel;
    }

    ApplyBoundaryConditions();
    CalculateParameters();
    //if (m_config.canonical) 
    //  RenormalizeVelocities();
    t += dt;
}


float*
MDSystem::getArray(int type)
{
    //assert(m_bInitialized);

    float* hdata = 0;
    float* ddata = 0;

    unsigned int pbo = 0;

    switch (type)
    {
    default:
    case 0:
        hdata = h_Pos;
        ddata = d_Pos;
    //if (usepbo)
    //    pbo = m_vbo[currentRead];
        break;
    case 1:
        hdata = h_Vel;
        ddata = d_Vel;
        break;
    }

    copyArrayFromDevice(hdata, ddata, pbo, m_config.N);
    return hdata;
}

void
MDSystem::setArray(int type, const float* data)
{
    //assert(m_bInitialized);

    switch (type)
    {
    default:
    case 0:
        copyArrayToDevice(d_Pos, data, m_config.N);
        break;
    case 1:
        copyArrayToDevice(d_Vel, data, m_config.N);
        break;
    }
}


SplineFunction MDSystem::RDF(double rmax, double shag)
{
    const int &N = m_config.N;
    {
      double n0 = N / L / L / L;
      std::vector<double> rdf_x, rdf_y;
      for (int ir = 0; ir < NdNdr2.size(); ++ir) {
        double r2 = (ir + 0.5) * rdf_dr2;
        double r = sqrt(r2);
        rdf_x.push_back(r);
        rdf_y.push_back(NdNdr2[ir] / rdf_dr2 / 2. / PI / r / n0 / static_cast<double>(N));
      }

      SplineFunction ret(rdf_x, rdf_y);
      return ret;
    }
}

void MDSystem::initvelo(double vmax, double step)
{
    const int &N = m_config.N;
    int dens[100000];
    int maxind = (int)(vmax/ step) + 1;
    int tind = 0;
    memset(dens, 0, sizeof(dens));
    for(int i=0;i<4*N;i+=4)
    {
        tind = (int)(sqrt(h_Vel[i]*h_Vel[i]+h_Vel[i+1]*h_Vel[i+1]+h_Vel[i+2]*h_Vel[i+2]) / step);
        if (tind<maxind) dens[tind]++;
    }
    curvelo = SplineFunction();
    for(int i=0;i<maxind;++i)
    {
        curvelo.add_val(step *(0.5 + i), dens[i] / step / N);
    }
    veloIters = 1;
}

SplineFunction MDSystem::getvelo()
{
    return curvelo;
}

void MDSystem::updatevelo()
{
    const int &N = m_config.N;
    int dens[100000];
    int maxind = curvelo.vals.size();//(int)(vmax/shag) + 1;
    int tind = 0;
    double shag = curvelo.vals[1].first - curvelo.vals[0].first;
    memset(dens, 0, sizeof(dens));
    for(int i=0;i<4*N;i+=4)
    {
        tind = (int)(sqrt(h_Vel[i]*h_Vel[i]+h_Vel[i+1]*h_Vel[i+1]+h_Vel[i+2]*h_Vel[i+2]) / shag );
        if (tind<maxind) dens[tind]++;
    }
    for(int i=0;i<maxind;++i)
    {
        curvelo.vals[i].second = (curvelo.vals[i].second*veloIters + dens[i] / shag / N) / (veloIters + 1);
    }
    veloIters++;
}

double MDSystem::Maxwell(double v)
{
    return 4*PI*pow(/*m_config.mass*/1./2./PI/T, 3./2.)*v*v*exp(-/*m_config.mass**/v*v/2./T);
}

void MDSystem::resetAveraging()
{
  av_U_tot = av_T_tot = av_p_tot = 0.;
  av_iters = 0;
}

// Fluctuations in a subvolume compared to the binomial distribution
void MDSystem::Fluctuations(double fraction)
{
    const int &N = m_config.N;
    double tsz = L * pow(fraction, 1./3.);
    long long tN = 0;
    for(int i=0;i<4*N;i+=4)
    {
        if (h_Pos[i]>=0.5*L-0.5*tsz && h_Pos[i]<=0.5*L+0.5*tsz
                && h_Pos[i+1]>=0.5*L-0.5*tsz && h_Pos[i+1]<=0.5*L+0.5*tsz
                && h_Pos[i+2]>=0.5*L-0.5*tsz && h_Pos[i+2]<=0.5*L+0.5*tsz)
            tN++;
    }
    momN  += tN;
    momN2 += tN*tN;
    momN3 += tN*tN*tN;
    momN4 += tN*tN*tN*tN;
    momiters++;

    std::cout << "Iteration: " << momiters << "\t" << "<N> = " << momN / momiters
            << "\t" << "w[N] = " << (momN2/momiters - (momN/momiters)*(momN/momiters)) / (momN/momiters)
               << "\t" << "EV: " << (1. - 2.*PI/3. * N / L / L / L) * (1. - 2.*PI/3. * N / L / L / L)
                  << "\t" << "Binom: " << 1. - fraction << "\n";
}

int MDSystem::fast_round(float x)
{
	if (x > 0)
		return static_cast<int>(x + 0.5f);
	else
		return static_cast<int>(x - 0.5f);
	//return lround(x);
}

void MDSystem::setHardwareMode(bool useCUDA)
{
  m_config.useCUDA = useCUDA;
  ReallocateMemory();
}
