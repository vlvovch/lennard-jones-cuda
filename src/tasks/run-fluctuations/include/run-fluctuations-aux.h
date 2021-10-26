#ifndef RUNFLUCTUATIONSAUX_H
#define RUNFLUCTUATIONSAUX_H

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


struct RunFluctuationsParameters {
  // Output file names prefix
  std::string output_prefix;

  // Numeric parameters
  std::map<std::string, double> parameters;

  // Default constructor with the default values
  RunFluctuationsParameters(const std::string& input_file = "") :
    output_prefix("run"),
    parameters({
      {"N",           400},        // the number of particles
      {"T*",          1.4},        // the temperature, ignored if the microcanonical ensemble is used
      //{"u*",          1.4},        // the energy per particle, ignored if the canonical ensemble is used
      {"u*",          1.708},        // the energy per particle, ignored if the canonical ensemble is used
      {"rho*",        0.05},       // the particle density
      {"teq",         50.},        // the equilibration time time, if negative, the calculation goes on indefinitely until stopped
      {"tfin",        1000.},      // the maximum time, 
      {"dt*",         0.004},      // the integration time step
      {"canonical",   0},          // the ensemble: 0 - microcanonical, 1 - canonical
      {"subvolume_spacing", 0.05}, // the spacing in the values of the subvolume fractions
      {"useCUDA",     1}           // if available, use CUDA GPU
      })
  {
    if (input_file != "")
      ReadParametersFromFile(input_file);
  }

  void ReadParametersFromFile(const std::string& filename) {
    std::ifstream fin(filename);
    if (fin.is_open()) {
      std::cout << "Reading input parameters from file " << filename << "\n";
      std::string var;
      double val;

      while (fin >> var) {
        if (var.size() == 0 || var[0] == '#') {
          fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
          continue;
        }

        std::cout << "Reading input parameter " << var << " = ";

        if (var == "output_prefix") {
          fin >> output_prefix;
          std::cout << output_prefix << std::endl;
          continue;
        }

        fin >> val;
        std::cout << val << std::endl;

        parameters[var] = val;
      }
    }
    else {
      std::cout << "Cannot open parameters file!" << "\n";
    }
    std::cout.flush();
  }

  void OutputParameters() {
    const int tabsize = 25;

    std::cout << "Lennard-Jones Molecular Dynamics fluctuations run parameter list:" << "\n";

    for (auto& el : parameters) {
      std::cout << std::setw(tabsize) << el.first << " = " << el.second << "\n";
    }

    std::cout << std::setw(tabsize) << "output_prefix" << " = " << output_prefix << "\n";

    std::cout << std::endl;
  }

  std::string GetFullPrefix() {
    std::string ret = output_prefix;

    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::stringstream ss;
    ss << std::put_time(&tm, "%d-%m-%Y-T%H-%M-%S");
    ss << ".N" << parameters["N"];
    if (std::lround(parameters["canonical"]))
      ss << ".Tst" << parameters["T*"];
    else
      ss << ".ust" << parameters["u*"];

    ss << ".rhost" << parameters["rho*"];

    ret += "." + ss.str();

    return ret;
  }

};

namespace RunFluctuationsFunctions {

  // Calculate the number of particles in a subsystem at a given moment
  // type: 0 - |x| < xcut, 1 - |y| < ycut, 2 - |z| < zcut, 3 - cube around the center
  // fraction: fraction of the total volume
  static int GetNSubsystem(const MDSystem& syst, double fraction = 0.5, int type = 3) {
    double L = syst.L;
    int N = syst.m_config.N;

    double tL = L * pow(fraction, 1. / 3.);

    int ret = 0;

    for (int i = 0; i < 4 * N; i += 4) {
      double x = syst.h_Pos[i];
      double y = syst.h_Pos[i + 1];
      double z = syst.h_Pos[i + 2];
      if (type == 3) {
        double dx = (x - 0.5 * L);
        double dy = (y - 0.5 * L);
        double dz = (z - 0.5 * L);

        if (abs(dx) < 0.5 * tL && abs(dy) < 0.5 * tL && abs(dz) < 0.5 * tL)
          ret++;
      }
      else {
        double coord = x;
        if (type == 1)
          coord = y;
        else if (type == 2)
          coord = z;
        if (coord < fraction * L)
          ret++;
      }
    }

    return ret;
  }

  // Calculate the number of particles in a momentum space subsystem at a given time moment
  // type: 0 - |vx| < vxcut, 1 - |vy| < vycut, 2 - |vz| < vzcut
  // fraction: fraction of the total volume
  static int GetNsubVz(const MDSystem& syst, double vcut = 0.5, int type = 0) {
    int N = syst.m_config.N;

    int ret = 0;

    for (int i = 0; i < 4 * N; i += 4) {
      double vx = syst.h_Vel[i];
      double vy = syst.h_Vel[i + 1];
      double vz = syst.h_Vel[i + 2];

      double v = vx;
      if (type == 1)
        v = vy;
      else if (type == 2)
        v = vz;

      if (abs(v) < vcut)
        ret++;
    }

    return ret;
  }

  // Calculate the number of particles in subsystems at a given time moment
  // type: 0 - |x| < xcut, 1 - |y| < ycut, 2 - |z| < zcut, 3 - cube around the center
  // alpha_step: step in the fractions of the total volume, from 0 to 1
  static std::vector<int> GetNSubsystemBatch(const MDSystem& syst, double alpha_step = 0.05, int type = 3) {
    double L = syst.L;
    int N = syst.m_config.N;

    std::vector<double> tLs;
    std::vector<int> cnts;
    
    for (double alpha = alpha_step; alpha < 1. - 1.e-9; alpha += alpha_step) {
      if (type == 3) {
        tLs.push_back(L * pow(alpha, 1. / 3.));
      }
      cnts.push_back(0);
    }

    for (int i = 0; i < 4 * N; i += 4) {
      double x = syst.h_Pos[i];
      double y = syst.h_Pos[i + 1];
      double z = syst.h_Pos[i + 2];
      if (type == 3) {
        double dx = (x - 0.5 * L);
        double dy = (y - 0.5 * L);
        double dz = (z - 0.5 * L);

        //if (abs(dx) < 0.5 * tL && abs(dy) < 0.5 * tL && abs(dz) < 0.5 * tL)
        //  ret++;
        int tindx = std::distance(tLs.begin(), upper_bound(tLs.begin(), tLs.end(), abs(dx) / 0.5));
        int tindy = std::distance(tLs.begin(), upper_bound(tLs.begin(), tLs.end(), abs(dy) / 0.5));
        int tindz = std::distance(tLs.begin(), upper_bound(tLs.begin(), tLs.end(), abs(dz) / 0.5));
        int tind = std::max({ tindx,tindy,tindz });
        if (tind < cnts.size())
          cnts[tind]++;
      }
      else {
        double coord = x;
        if (type == 1)
          coord = y;
        else if (type == 2)
          coord = z;

        coord /= L;
        int tind = (int)(coord / alpha_step);
        if (tind < cnts.size())
          cnts[tind]++;
      }
    }

    for (int i = 1; i < cnts.size(); ++i)
      cnts[i] += cnts[i - 1];

    return cnts;
  }

  // Calculate the number of particles in a momentum space subsystem at a given time moment
  // type: 0 - |vx| < vxcut, 1 - |vy| < vycut, 2 - |vz| < vzcut
  // alpha_step: step in the fractions of the total volume, from 0 to 1
  static std::vector<int> GetNsubVzBatch(const MDSystem& syst, double vcut_max, double alpha_step = 0.05, int type = 0) {
    int N = syst.m_config.N;

    std::vector<int> cnts;

    for (double alpha = alpha_step; alpha < 1. + 1.e-9; alpha += alpha_step) {
      cnts.push_back(0);
    }

    int ret = 0;

    for (int i = 0; i < 4 * N; i += 4) {
      double vx = syst.h_Vel[i];
      double vy = syst.h_Vel[i + 1];
      double vz = syst.h_Vel[i + 2];

      double v = vx;
      if (type == 1)
        v = vy;
      else if (type == 2)
        v = vz;

      v /= vcut_max;
      int tind = (int)(abs(v) / alpha_step);
      if (tind < cnts.size())
        cnts[tind]++;
    }

    for (int i = 1; i < cnts.size(); ++i)
      cnts[i] += cnts[i - 1];

    return cnts;
  }

  // Get the average velocities
  // Returns a vector {<vx>, <vy>, <vz>}
  static std::vector<double> GetAvVel(const MDSystem& syst) {
    std::vector<double> ret(3, 0.);
    for (int i = 0; i < 4 * syst.m_config.N; i += 4) {
      ret[0] += syst.h_Vel[i];
      ret[1] += syst.h_Vel[i + 1];
      ret[2] += syst.h_Vel[i + 2];
    }
    for (int i = 0; i < 3; ++i)
      ret[i] /= syst.m_config.N;
    return ret;
  }


  // Radial distribution function time average
  struct RDF_average {
    double max_r;
    double dr;
    int iters_rdf;
    SplineFunction Gr;

    RDF_average(double maxr = 5.0, double ddr = 0.03) : max_r(maxr), dr(ddr), iters_rdf(0) {}

    void AddTimeStep(MDSystem& syst) {
      SplineFunction crdf = syst.RDF(max_r, dr);
      if (iters_rdf == 0)
        Gr = crdf;
      else
      {
        for (int i = 0; i < crdf.vals.size(); ++i)
          Gr.vals[i].second = (Gr.vals[i].second * iters_rdf + crdf.vals[i].second) / (iters_rdf + 1);
      }
      iters_rdf++;
    }

    void PrintToFile(const std::string& filename) {
      std::ofstream fout(filename);
      if (fout.is_open()) {
        fout << std::setw(15) << "r*" << " ";
        fout << std::setw(15) << "G(r)" << " ";
      }
      fout << std::endl;

      for (int i = 0; i < Gr.vals.size(); ++i) {
        fout << std::setw(15) << Gr.vals[i].first << " ";
        fout << std::setw(15) << Gr.vals[i].second << " ";
        fout << std::endl;
      }

      fout.close();
    }
  };

  // Time average for coordinate space fluctuations
  struct CoordFlucsAverage {
    int type;
    double alpha_step;
    std::vector<double> alphas, totsN, totsN2;
    int iters;

    CoordFlucsAverage(int type_ = 0, double alphastep = 0.05) :
      type(type_),
      alpha_step(alphastep),
      iters(0) 
    {
      totsN = totsN2 = std::vector<double>(0);
      for (double alpha = alphastep; alpha <= 1. - 1.e-9; alpha += alpha_step) {
        totsN.push_back(0.);
        totsN2.push_back(0.);
        alphas.push_back(alpha);
      }
    }

    void AddTimeStep(const MDSystem& syst) {
      std::vector<int> cnts = GetNSubsystemBatch(syst, alpha_step, type);
      for (int ia = 0; ia < alphas.size(); ++ia) {
        double alpha = alphas[ia];
        //int Nsub = GetNSubsystem(syst, alpha, type);
        int Nsub = cnts[ia];
        totsN[ia]  += static_cast<double>(Nsub);
        totsN2[ia] += static_cast<double>(Nsub) * static_cast<double>(Nsub);

        /*if (ia == alphas.size() - 1) {
          std::cout << Nsub << " " << cnts[ia] << " ; ";
        }*/
      }
      iters++;
    }

    void PrintToFile(const std::string& filename) {
      std::ofstream fout(filename);
      if (fout.is_open()) {
        fout << std::setw(15) << "alpha" << " ";
        fout << std::setw(15) << "<N>" << " ";
        fout << std::setw(15) << "w[N]" << " ";
        fout << std::setw(15) << "w[N]/(1-alpha)" << " ";
      }
      fout << std::endl;

      for (int ia = 0; ia < alphas.size(); ++ia) {
        double alpha = alphas[ia];

        fout << std::setw(15) << alpha << " ";

        double Nav  = totsN[ia] / iters;
        double N2av = totsN2[ia] / iters;
        double w = (N2av - Nav * Nav) / Nav;

        fout << std::setw(15) << Nav << " ";
        fout << std::setw(15) << w << " ";
        fout << std::setw(15) << w / (1. - alpha) << " ";
        fout << std::endl;
      }

      fout.close();
    }
  };


  struct MomentumFlucsAverage {
    int type;
    double step;
    const double Vfactor = 2.0;
    std::vector<double> vcuts, totsN, totsN2;
    int iters;
    int totN;

    MomentumFlucsAverage(const MDSystem& syst, int type_ = 0, double alphastep = 0.05) :
      type(type_),
      step(alphastep),
      iters(0)
    {
      totsN = totsN2 = std::vector<double>(0);
      double Tst = syst.m_config.T0;
      totN = syst.m_config.N;
      for (double alpha = alphastep; alpha <= 1. + 1.e-9; alpha += step) {
        double vcut = sqrt(Tst) * alpha * Vfactor;
        totsN.push_back(0.);
        totsN2.push_back(0.);
        vcuts.push_back(vcut);
      }
    }

    void AddTimeStep(const MDSystem& syst) {
      std::vector<int> cnts = GetNsubVzBatch(syst, sqrt(syst.m_config.T0) * Vfactor, step, type);
      for (int ia = 0; ia < vcuts.size(); ++ia) {
        double vcut = vcuts[ia];
        //int Nsub = GetNsubVz(syst, vcut, type);
        int Nsub = cnts[ia];
        totsN[ia] += static_cast<double>(Nsub);
        totsN2[ia] += static_cast<double>(Nsub) * static_cast<double>(Nsub);

        //if (ia == vcuts.size() - 1) {
        //  std::cout << Nsub << " " << cnts[ia] << " ; ";
        //}
      }
      iters++;
    }

    void PrintToFile(const std::string& filename) {
      std::ofstream fout(filename);
      if (fout.is_open()) {
        fout << std::setw(15) << "vcut" << " ";
        fout << std::setw(15) << "<N>" << " ";
        fout << std::setw(15) << "alpha" << " ";
        fout << std::setw(15) << "w[N]" << " ";
        fout << std::setw(15) << "w[N]/(1-alpha)" << " ";
      }
      fout << std::endl;

      for (int ia = 0; ia < vcuts.size(); ++ia) {
        double vcut = vcuts[ia];

        fout << std::setw(15) << vcut << " ";

        double Nav = totsN[ia] / iters;
        double N2av = totsN2[ia] / iters;
        double w = (N2av - Nav * Nav) / Nav;

        fout << std::setw(15) << Nav << " ";
        double alpha = Nav / totN;
        fout << std::setw(15) << alpha << " ";
        fout << std::setw(15) << w << " ";
        fout << std::setw(15) << w / (1. - alpha) << " ";
        fout << std::endl;
      }

      fout.close();
    }
  };
}
#endif