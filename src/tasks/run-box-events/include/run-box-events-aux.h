#ifndef RUNBOXEVENTSAUX_H
#define RUNBOXEVENTSAUX_H

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
#include <ostream>


#ifdef __linux__
#include <sys/stat.h>
#endif

#include "MDSystem.h"
#include "NumberStatistics.h"


struct RunBoxEventsParameters {
  // Output file names prefix
  std::string output_prefix;

  // Numeric parameters
  std::map<std::string, double> parameters;

  // Default constructor with the default values
  RunBoxEventsParameters(const std::string& input_file = "") :
    output_prefix("run"),
    parameters({
      {"nevents",    1000},        // the number of events
      {"N",           400},        // the number of particles
      {"T*",          1.4},        // the temperature, ignored if the microcanonical ensemble is used
      {"u*",          -0.0458},    // the energy per particle, ignored if the canonical ensemble is used
      {"rho*",        0.30},       // the particle density
      {"dt_out",      2.},         // the output interval time
      {"tfin",        50.},        // the maximum time per single simulation
      {"dt*",         0.004},      // the integration time step
      {"canonical",   0},          // the ensemble: 0 - microcanonical, 1 - canonical
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

    std::cout << "Lennard-Jones Molecular Dynamics box simulation event generator parameter list:" << "\n";

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

namespace RunBoxEventsFunctions {
  const int tabsize = 15;

  void OutputEvent(const MDSystem& syst, std::ostream& out, int nevent = 1) {
    out << "Event # " << nevent << "\n";
    out << "   N = " << syst.m_config.N << "\n";
    out << "  u* = " << syst.U / syst.m_config.N << "\n";
    out << "  T* = " << syst.T << "\n";
    out << "rho* = " << syst.m_config.N / pow(syst.L,3) << "\n";
    out << "  t* = " << syst.t << "\n";
    out << std::setw(tabsize) << "x" << " ";
    out << std::setw(tabsize) << "y" << " ";
    out << std::setw(tabsize) << "z" << " ";
    out << std::setw(tabsize) << "vx" << " ";
    out << std::setw(tabsize) << "vy" << " ";
    out << std::setw(tabsize) << "vz" << " ";
    out << std::endl;
    for(int iN = 0; iN < 4*syst.m_config.N; iN += 4) {
      out << std::setw(tabsize) << syst.h_Pos[iN] << " ";
      out << std::setw(tabsize) << syst.h_Pos[iN + 1] << " ";
      out << std::setw(tabsize) << syst.h_Pos[iN + 2] << " ";
      out << std::setw(tabsize) << syst.h_Vel[iN] << " ";
      out << std::setw(tabsize) << syst.h_Vel[iN + 1] << " ";
      out << std::setw(tabsize) << syst.h_Vel[iN + 2] << " ";
      out << std::endl;
    }
    out << std::endl;
  }

  void OutputEvent(const MDSystem& syst, const std::string& fileout, int nevent = 1) {
    std::ofstream fout(fileout, std::ios::app);
    OutputEvent(syst, fout, nevent);
    fout.close();
  }
} // namespace RunBoxEventsFunctions

#endif