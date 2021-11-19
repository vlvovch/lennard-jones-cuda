#ifndef RUNISOTHERMAUX_H
#define RUNISOTHERMAUX_H

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

#include "../../auxiliary/time-average-aux.h"


struct RunIsothermParameters {
  // Output file names prefix
  std::string output_prefix;

  // Numeric parameters
  std::map<std::string, double> parameters;

  // Default constructor with the default values
  RunIsothermParameters(const std::string& input_file = "") :
    output_prefix("isotherm.run"),
    parameters({
      {"N",           400},        // the number of particles
      {"T*",          1.4},        // the temperature of the isotherm
      {"rho*_min",    0.60},       // the minimum particle density
      {"rho*_max",    0.60},       // the maximum particle density
      {"drho*",       0.01},       // the step in particle density
      {"teq",         10.},        // the equilibration time time
      {"tfin",        5000.},      // the maximum time per single simulation
      {"dt*",         0.004},      // the integration time step
      //{"canonical",   1},          // the ensemble: 0 - microcanonical, 1 - canonical
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

    std::cout << "Lennard-Jones Molecular Dynamics isotherm run parameter list:" << "\n";

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
    //if (std::lround(parameters["canonical"]))
      ss << ".Tst" << parameters["T*"];
    //else
    //  ss << ".ust" << parameters["u*"];

    //ss << ".rhost" << parameters["rho*"];

    ret += "." + ss.str();

    return ret;
  }

};

namespace RunIsothermFunctions {


}

#endif