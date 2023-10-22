#include "NumberStatistics.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace SampleMoments;

const int tabsize = 15;

std::string filenames_file = "input.txt";
double dalpha = 0.1;
int Nalphas = 10;
double L = 1., dL = 0.1;
double t = 50.;

double alpha_i(int i) {
  return dalpha * i;
}

int main(int argc, char *argv[]) {

  if (argc > 1) {
    filenames_file = string(argv[1]);
  }

  if (argc > 2) {
    dalpha = atof(argv[2]);
  }

  Nalphas = round(1. / dalpha) + 1;
  if ((Nalphas-1) * dalpha > 1.)
    Nalphas--;

  cout << Nalphas << endl;

  vector<double> times;
  vector<string> filenames;
  // Read the filenames for each time moment
  {
    ifstream fin(filenames_file);
    if (fin.is_open()) {
      double t;
      string file;
      while (fin >> t >> file) {
        times.push_back(t);
        filenames.push_back(file);
      }
      fin.close();
    } else {
      cerr << "Error: could not open file " << filenames_file << endl;
      return 1;
    }
  }

  ofstream fout( filenames_file + ".fluctuations-time-dep.dat");
  if (!fout.is_open()) {
    cerr << "Error: could not open output file " << filenames_file << ".fluctuations-time-dep.dat" << endl;
    return 1;
  }

  fout << Nalphas << " # Number of subvolumes" << endl;
  fout << endl;

  cout   << setw(tabsize) << "t" << " "
         << setw(tabsize) << "nevents" << " "
         << setw(tabsize) << "alpha" << " "
         << setw(tabsize) << "mean" << " "
         << setw(tabsize) << "error" << " "
         << setw(tabsize) << "wtil" << " "
         << setw(tabsize) << "error" << " ";
//         << setw(tabsize) << "Stil" << " "
//         << setw(tabsize) << "error" << " ";
  cout << endl;

  for(int ifile = 0; ifile < filenames.size(); ifile++) {
    const string &filename = filenames[ifile];
    double t = times[ifile];

    fout << t << " # Time" << endl;

    vector<NumberStatistics> stats(Nalphas);

    long long nevents = 0;

    double rho = 0.1;
    int N = 400;

    ifstream fin(filename);
    string in1, in2;
    int in_int;
    double in_double;

    if (fin.is_open()) {
      // Read the events
      while (fin >> in1 >> in2 >> in_int) {
        // Read single event
        assert(in1 == "Event" && in2 == "#");
        nevents++;

//        if (nevents % 100 == 0)
//          cout << "Reading event #" << nevents << endl;

        fin >> in1 >> in2 >> in_int;
        assert(in1 == "N" && in2 == "=");
        if (nevents == 1) {
          N = in_int;
          //cout << "  N = " << N << endl;
        } else {
          assert(N == in_double);
        }
        fin >> in1 >> in2 >> in_double;
        assert(in1 == "u*" && in2 == "=");
        fin >> in1 >> in2 >> in_double;
        assert(in1 == "T*" && in2 == "=");
        fin >> in1 >> in2 >> in_double;
        assert(in1 == "rho*" && in2 == "=");
        if (nevents == 1) {
          rho = in_double;
          //cout << "rho* = " << rho << endl;
          L = pow(N / rho, 1. / 3.);
          //cout << "  L* = " << L << endl;
        } else {
          assert(rho == in_double);
        }
        fin >> in1 >> in2 >> in_double;
        assert(in1 == "t*" && in2 == "=");
        if (nevents == 1) {
          t = in_double;
          //cout << "  t* = " << t << endl;
        } else {
          assert(t == in_double);
        }

        // Header line
        getline(fin, in1);
        getline(fin, in1);

        // Initialize counters
        vector<int> cnts(Nalphas);

        // Process particles
        double x, y, z, vx, vy, vz;
        for (int iN = 0; iN < N; iN++) {
          fin >> x >> y >> z >> vx >> vy >> vz;

          // Get the bin and count particle
          int indz = static_cast<int>((z / L) / dalpha) + 1;
          if (indz < Nalphas)
            cnts[indz]++;
        }

        // Compute the prefix sums
        for (int i = 1; i < Nalphas; i++)
          cnts[i] += cnts[i - 1];

        // Add the counts to the statistics
        for (int i = 0; i < Nalphas; i++)
          stats[i].AddObservation(cnts[i]);
      }
    } else {
      cout << "Cannot open file " << filename << endl;
      return 1;
    }

//    // Print the statistics
//    cout << "Statistics for " << nevents << " events" << endl;
//    cout << setw(tabsize) << "alpha" << " "
//         << setw(tabsize) << "mean" << " "
//         << setw(tabsize) << "error" << " "
//         << setw(tabsize) << "wtil" << " "
//         << setw(tabsize) << "error" << " "
//         << setw(tabsize) << "Stil" << " "
//         << setw(tabsize) << "error" << " ";
//    cout << endl;
//
//    for (int i = 0; i < Nalphas; i++) {
//      double alpha = alpha_i(i);
//      cout << setw(tabsize) << alpha << " "
//           << setw(tabsize) << stats[i].GetMean() << " "
//           << setw(tabsize) << stats[i].GetMeanError() << " "
//           << setw(tabsize) << stats[i].GetScaledVariance() / (1. - alpha) << " "
//           << setw(tabsize) << stats[i].GetScaledVarianceError() / (1. - alpha) << " "
//           << setw(tabsize) << stats[i].GetSkewness() / abs(1. - 2. * alpha) << " "
//           << setw(tabsize) << stats[i].GetSkewnessError() / abs(1. - 2. * alpha) << " ";
//      cout << endl;
//    }

    {
      int i = Nalphas / 2;
      double alpha = alpha_i(i);
      cout << setw(tabsize) << t << " "
           << setw(tabsize) << nevents << " "
           << setw(tabsize) << alpha << " "
           << setw(tabsize) << stats[i].GetMean() << " "
           << setw(tabsize) << stats[i].GetMeanError() << " "
           << setw(tabsize) << stats[i].GetScaledVariance() / (1. - alpha) << " "
           << setw(tabsize) << stats[i].GetScaledVarianceError() / (1. - alpha) << " ";
//           << setw(tabsize) << stats[i].GetSkewness() / abs(1. - 2. * alpha) << " "
//           << setw(tabsize) << stats[i].GetSkewnessError() / abs(1. - 2. * alpha) << " ";
      cout << endl;
    }

    fout << nevents << " # Number of events" << endl;
    fout << setw(tabsize) << "alpha" << " "
         << setw(tabsize) << "mean" << " "
         << setw(tabsize) << "error" << " "
         << setw(tabsize) << "wtil" << " "
         << setw(tabsize) << "error" << " "
         << setw(tabsize) << "Stil" << " "
         << setw(tabsize) << "error" << " ";
    fout << endl;

    for (int i = 0; i < Nalphas; i++) {
      double alpha = alpha_i(i);
      fout << setw(tabsize) << alpha << " "
           << setw(tabsize) << stats[i].GetMean() << " "
           << setw(tabsize) << stats[i].GetMeanError() << " "
           << setw(tabsize) << stats[i].GetScaledVariance() / (1. - alpha) << " "
           << setw(tabsize) << stats[i].GetScaledVarianceError() / (1. - alpha) << " "
           << setw(tabsize) << stats[i].GetSkewness() / abs(1. - 2. * alpha) << " "
           << setw(tabsize) << stats[i].GetSkewnessError() / abs(1. - 2. * alpha) << " ";
      fout << endl;
    }
    fout << endl;

  }

  fout.close();

  return 0;
}
 