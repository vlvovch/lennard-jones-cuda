#ifndef mdsystem_h
#define mdsystem_h
#include "splinefunction.h"
#include "MersenneTwister/MersenneTwister.h"

class MDSystem
{
public:

  struct MDSystemConfiguration {
    // Number of particles
    int N;

    // Initial/CE temperature
    double T0;

    // Density of particles
    double rho;

    // Canonical ensemble instead of microcanonical
    bool canonical;


    // Boundary conditions. 0 - periodic, 1 - hard wall
    int boundaryConditions;

    // Use CUDA GPU (if available)
    bool useCUDA;
    int CUDABlockSize;

    // Default values
    MDSystemConfiguration() {
      N = 128;
      T0 = 1.500;
      rho = 0.200;
      canonical = false;
      boundaryConditions = 0;
      useCUDA = false;
      CUDABlockSize = 256;
    }
  };

  MDSystemConfiguration m_config;

  class MaxwellDistributionGenerator
  {
    double probmax;
    MTRand rangen;
  public:
    MaxwellDistributionGenerator() {
      probmax = 1.2;
    }

    double prob(double alpha);

    double sampleAlphaV();
  };

  MaxwellDistributionGenerator m_MaxwellGenerator;

  // Dynamical parameters:
  // total energy, instant temperature, kinetic energy, potential energy, pressure
  double U, T, K, V, P;

  // CUDA initialized
  bool CUDAInit;

  // Vector with coordinats, momenta, and computed forces
  float* h_Pos;
  float* h_Vel;
  float* h_Force;
  float* d_Pos;
  float* d_Vel;
  float* d_Force;

  // Current time moment
  double t;

  // Side length
  double L;

  // Velocity distribution
  SplineFunction curvelo;
  int veloIters;

  // Radial distribution function
  std::vector<int> NdNdr2;
  float rdf_dr2;

  // Moments for fluctuations
  double momN, momN2, momN3, momN4;
  int momiters;

  // Time-averaged quantities
  double av_U_tot, av_T_tot, av_p_tot;
  int av_iters;

  // Random number generator
  MTRand rangen;
public:
  MDSystem(const MDSystem::MDSystemConfiguration& config = MDSystem::MDSystemConfiguration());
  
  ~MDSystem(void);

  void Reinitialize(const MDSystem::MDSystemConfiguration& config);

  void ReallocateMemory();

  void SampleInitialConditions();

  void CorrectTotalMomentum();

  void CalculateForces();

  void CalculateForces(float* frc);

  double KineticTemperature(float* vel);

  double CalculateXi(float* frc, float* vel); // For constant Tkin simulations

  void CalculateParameters();

  void RenormalizeVelocities(bool RecalculateTkin = false);

  void ApplyBoundaryConditions();

  double getTime() { return t; }

  void Integrate(double dt);

  // Get array from GPU memory
  float* getArray(int type);

  // Put array into GPU memory
  void setArray(int type, const float* data);

  // Velocity distribution
  void initvelo(double vmax = 12., double shag = 0.12);
  void updatevelo();//(double vmax = 2., double shag = 0.02);
  SplineFunction getvelo();
  double Maxwell(double v);

  // Reset averaging
  void resetAveraging();

  // Radial distribution function
  SplineFunction RDF(double rmax = 5., double shag = 0.1);

  // Fluctuations (currently not used)
  void Fluctuations(double fraction);

  // Rounding for the periodic boundary conditions
  int fast_round(float x);

  void setHardwareMode(bool useCUDA);
  void setCanonical(bool canonical) { m_config.canonical = canonical; }
  void setPeriodicBoundaryCondition(bool periodic) { m_config.boundaryConditions = !periodic; }
};

#endif
