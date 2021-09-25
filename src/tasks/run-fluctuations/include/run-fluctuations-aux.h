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
      {"u*",          1.4},        // the energy per particle, ignored if the canonical ensemble is used
      {"rho*",        0.05},       // the particle density
      {"teq",         50.},        // the equilibration time time, if negative, the calculation goes on indefinitely until stopped
      {"tfin",        1000.},      // the maximum time, 
      {"dt*",         0.004},      // the integration time step
      {"canonical",   1},          // the ensemble: 0 - microcanonical, 1 - canonical
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

  // Calculate the number of particles in a momentum space subsystem at a given moment
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
};
#endif