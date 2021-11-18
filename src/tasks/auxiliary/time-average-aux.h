#ifndef TIMEAVERAGEAUX_H
#define TIMEAVERAGEAUX_H

#include <fstream>
#include <string>
#include <functional>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <map>
#include <ctime>


#ifdef __linux__
#include <sys/stat.h>
#endif

#include "MDSystem.h"
#include "NumberStatistics.h"

// Calculate time average including stat. error. estimation
// via statistical inefficiency
// Assumes exponential decay of the correlator, 
// see pp. 194-195 in Allen, Tildesley, "Computer Simulation of Liquids"
struct TimeAverage
{
  SampleMoments::NumberStatistics stats;
  double lastN;
  double corrNtot;
  int iters;
  TimeAverage() {
    iters = 0;
    lastN = corrNtot = 0.;
  }

  void AddObservation(double val) {
    stats.AddObservation(val);

    if (iters > 0) {
      corrNtot += lastN * val;
    }

    lastN = val;

    iters++;
  }

  double GetS()  {
    double Nav = stats.GetMean();
    double Nvar = stats.GetVariance();
    double Ncorrav = corrNtot / (iters - 1) - Nav * Nav;
    double s = 2. / log(Nvar / Ncorrav);
    if (s < 0.)
      s = 1.;
    return s;
  }

  double GetMean() {
    return stats.GetMean();
  }

  double GetMeanError() {
    return stats.GetMeanError() * sqrt(GetS());
  }
};

#endif