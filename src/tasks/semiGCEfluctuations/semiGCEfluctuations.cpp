#include "MDSystem.h"
#include "NumberStatistics.h"

using namespace std;

int GetNsub(const MDSystem& syst, double fraction = 0.5) {
  double L = syst.L;
  int N = syst.m_config.N;

  double tL = L * pow(fraction, 1. / 3.);

  int ret = 0;

  for (int i = 0; i < 4 * N; i += 4) {
    double dx = (syst.h_Pos[i] - 0.5 * L);
    double dy = (syst.h_Pos[i + 1] - 0.5 * L);
    double dz = (syst.h_Pos[i + 2] - 0.5 * L);

    if (abs(dx) < 0.5 * tL && abs(dy) < 0.5 * tL && abs(dz) < 0.5 * tL)
      ret++;
  }

  return ret;
}

int main(int argc, char *argv[])
{
  int nev = 10000;
  vector<double> fracs;
  for (double tfr = 0.05; tfr <= 1.; tfr += 0.05)
    fracs.push_back(tfr);

  vector<NumberStatistics> StatsFlucs(fracs.size());

  int N = 512;// 4096;
  double Tst = 2.0;
  double rhost = 0.050;
  Tst = 1.312;
  rhost = 0.316;

  double dt = 0.005;
  int iterspreeq = 10000;// / 2;
  int itersstep = 200;// / 2;

  MDSystem::MDSystemConfiguration config;
  config.N = N;
  config.T0 = Tst;
  config.rho = rhost; 
  config.useCUDA = true;

  MDSystem syst(config);

  syst.Reinitialize(config);


  // Equilibration phase
  syst.m_config.canonical = true;
  for (int i = 0; i < iterspreeq; ++i) {
    syst.Integrate(dt);
    //syst.RenormalizeVelocities();
  }


  // Production phase
  syst.m_config.canonical = false;
  for (int iN = 1; iN <= nev; iN++) {
    //syst.Reinitialize(config);

    //for (int i = 0; i < iterspreeq; ++i) {
    //  syst.Integrate(dt);
    //  syst.RenormalizeVelocities();
    //}

	//syst.RenormalizeVelocities();

    //syst.m_config.canonical = true;

    for (int i = 0; i < itersstep; ++i) {
      syst.Integrate(dt);
    }

    for (int i = 0; i < fracs.size(); ++i) {
      StatsFlucs[i].AddEvent(static_cast<double>(GetNsub(syst, fracs[i])));
    }

    //printf("%d ", iN);
    if (iN % 10 == 0) {
      for (int i = 0; i < fracs.size(); ++i) {
        StatsFlucs[i].CalculateCentralMoments();
        printf("%15d %10lf +- %-10lf %10lf +- %-10lf %10lf +- %-10lf %10lf +- %-10lf\n",
          iN,
          StatsFlucs[i].GetMean(),
          StatsFlucs[i].GetMeanError(),
          StatsFlucs[i].GetScaledVariance(),
          StatsFlucs[i].GetScaledVarianceError(), 
          StatsFlucs[i].GetScaledVariance() / (1. - fracs[i]),
          StatsFlucs[i].GetScaledVarianceError() / (1. - fracs[i]),
		      StatsFlucs[i].GetSkewness() / (1. - 2. * fracs[i]),
		      StatsFlucs[i].GetSkewnessError() / (1. - 2. * fracs[i])
          );
      }
      printf("\n");
    }
  }
  return 0;
}
 