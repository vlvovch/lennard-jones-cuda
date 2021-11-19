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


struct AnalyzeExpDecay
{
	SampleMoments::NumberStatistics mean;
	std::vector<SampleMoments::NumberStatistics> stats;
	std::vector<int> intervals;
	AnalyzeExpDecay(int dit = 100, int maxiters = 1000) {
		for (int diter = 0; diter <= maxiters; diter += dit) {
			intervals.push_back(diter);
		}
		stats.resize(intervals.size());
	}
	void AddObservation(const std::vector<double>& vals) {
		for (int ii = 0; ii < intervals.size(); ++ii) {
			if (vals.size() > intervals[ii]) {
				stats[ii].AddObservation(vals[vals.size() - 1] * vals[vals.size() - 1 - intervals[ii]]);
			}
		}
		mean.AddObservation(vals.back());
	}

	void AddSegment(const std::vector<double>& vals) {
		for (int ii = 0; ii < intervals.size(); ++ii) {
			if (vals.size() > intervals[ii]) {
				stats[ii].AddObservation(vals[0] * vals[intervals[ii]]);
			}
		}

		mean.AddObservation(vals[0]);
	}

	double GetTau(int interval) {
		if (interval <= 0. || interval >= intervals.size())
			return 0.;
		return intervals[interval] / log(mean.GetVariance() / (stats[interval].GetMean() - mean.GetMean() * mean.GetMean()));
	}

	double GetTauLocal(int interval) {
		if (interval <= 0. || interval >= intervals.size())
			return 0.;
		return (intervals[interval] - intervals[interval - 1]) / log((stats[interval - 1].GetMean() - mean.GetMean() * mean.GetMean()) / (stats[interval].GetMean() - mean.GetMean() * mean.GetMean()));
	}

	double GetTauSquare(int interval) {
		if (interval <= 0. || interval >= intervals.size())
			return 0.;
		return sqrt(intervals[interval] * intervals[interval] / log(mean.GetVariance() / (stats[interval].GetMean() - mean.GetMean() * mean.GetMean())));
	}

	double GetTauIntegral() {
		double ret = 0.;
		double var = mean.GetVariance();
		for (int ii = 0; ii < intervals.size() - 1; ++ii) {
			if (stats[ii].GetNumberOfObservations() < 1000 || stats[ii + 1].GetNumberOfObservations() < 1000)
				continue;

			double corr1 = stats[ii].GetMean() - mean.GetMean() * mean.GetMean();
			double corr2 = stats[ii + 1].GetMean() - mean.GetMean() * mean.GetMean();

			ret += 0.5 * (corr1 + corr2) * (intervals[ii + 1] - intervals[ii]);
		}
		return ret / var;
	}

	double GetTauIntegral(int ind, double tau = -1.) {
		if (tau < 0.)
			tau = GetTauLocal(ind);
		double ret = 0.;
		double var = mean.GetVariance();
		for (int ii = 0; ii < ind; ++ii) {
			//if (stats[ii].GetNumberOfObservations() < 1000 || stats[ii + 1].GetNumberOfObservations() < 1000)
			//	continue;

			double corr1 = stats[ii].GetMean() - mean.GetMean() * mean.GetMean();
			double corr2 = stats[ii + 1].GetMean() - mean.GetMean() * mean.GetMean();

			ret += 0.5 * (corr1 + corr2) * (intervals[ii + 1] - intervals[ii]);
		}
		ret /= var;
		ret += (stats[ind].GetMean() - mean.GetMean() * mean.GetMean()) / var * tau;
		return ret;
	}

	void PrintToFile(const std::string& filename, double dt, double tau, double factor = 1.) {
		std::ofstream fout(filename);

		fout << std::setw(15) << "t" << " ";
		fout << std::setw(15) << "Acorr" << " ";
		fout << std::setw(15) << "Acorr_err" << " ";
		fout << std::setw(15) << "Acorr_norm" << " ";
		fout << std::setw(15) << "Acorr_norm_err" << " ";
		fout << std::setw(15) << "exp(-t/tau)" << " ";
		//fout << setw(15) << "exp(-n*t/tau)" << " ";
		fout << std::setw(15) << "tau_extr" << " ";
		fout << std::setw(15) << "tau_extr_loc" << " ";
		fout << std::setw(15) << "tau_integral" << " ";
		fout << std::setw(15) << "coefficient" << " ";
		//fout << setw(15) << "tau_sq_extr" << " ";
		fout << std::endl;

		for (int ii = 0; ii < intervals.size(); ++ii) {
			if (stats[ii].GetNumberOfObservations() == 0)
				continue;
			double t = intervals[ii] * dt;
			fout << std::setw(15) << t << " ";
			fout << std::setw(15) << stats[ii].GetMean() - mean.GetMean() * mean.GetMean() << " ";
			fout << std::setw(15) << stats[ii].GetMeanError() << " ";
			fout << std::setw(15) << (stats[ii].GetMean() - mean.GetMean() * mean.GetMean()) / mean.GetVariance() << " ";
			fout << std::setw(15) << stats[ii].GetMeanError() / mean.GetVariance() << " ";
			fout << std::setw(15) << exp(-t / tau) << " ";
			//fout << setw(15) << exp(-intervals[ii] * t / tau) << " ";
			if (ii != 0)
				fout << std::setw(15) << dt * GetTau(ii) << " ";
			else
				fout << std::setw(15) << tau << " ";

			if (ii != 0)
				fout << std::setw(15) << dt * GetTauLocal(ii) << " ";
			else
				fout << std::setw(15) << tau << " ";

			if (ii != 0)
				fout << std::setw(15) << dt * GetTauIntegral(ii, GetTauLocal(ii)) << " ";
			else
				fout << std::setw(15) << dt * GetTauIntegral(ii, tau / dt) << " ";

			if (ii != 0)
				fout << std::setw(15) << dt * GetTauIntegral(ii, GetTauLocal(ii)) * mean.GetVariance() * factor << " ";
			else
				fout << std::setw(15) << dt * GetTauIntegral(ii, tau / dt) * mean.GetVariance() * factor << " ";

			//if (ii != 0)
			//	fout << setw(15) << dt * GetTauSquare(ii) << " ";
			//else
			//	fout << setw(15) << tau << " ";
			fout << std::endl;
		}

		fout << std::endl;
		//fout << "tau_int = " << GetTauIntegral() * dt << endl;
	}
};

#endif