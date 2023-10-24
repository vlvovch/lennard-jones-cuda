#include "include/run-box-events-aux.h"
#include "NumberStatistics.h"
#include "Utility.h"
#ifdef USE_CUDA_TOOLKIT
#include "cuda_runtime.h"
#endif

using namespace std;

const int tabsize = 15;

int main(int argc, char *argv[])
{
	RunBoxEventsParameters params;

	if (argc > 1) {
		cout << "Reading parameters from file " << string(argv[1]) << endl;
		params.ReadParametersFromFile(string(argv[1]));
	}
	params.output_prefix = params.GetFullPrefix();

	// check CUDA
	{
		int deviceCount = 0;
#ifdef USE_CUDA_TOOLKIT
		cudaGetDeviceCount(&deviceCount);
#endif

		if (deviceCount == 0 && lround(params.parameters["useCUDA"])) {
			std::cout << "Could not find a CUDA device! Calculations will be performed on CPU\n";
			params.parameters["useCUDA"] = 0;
		}
	}

	// Print all the run parameters
	params.OutputParameters();

	//cout << params.GetFullPrefix() << endl;

  int nevents = lround(params.parameters["nevents"]);
	int N = lround(params.parameters["N"]);
	double rhost = params.parameters["rho*"];
	double Tst = params.parameters["T*"];
	double ust = params.parameters["u*"];

	// Time parameters
	double dt   = params.parameters["dt*"];
	double tfin = params.parameters["tfin"];

	double dt_out  = params.parameters["dt_out"];

	MDSystem::MDSystemConfiguration config;
	config.N = N;
	config.T0 = Tst;
	config.rho = rhost; 
	config.useCUDA = lround(params.parameters["useCUDA"]);
	config.canonical = lround(params.parameters["canonical"]);

	MDSystem syst(config);

  // Do not compute RDF (~%20 speed-up)
  syst.m_computeRDF = false;
  syst.m_randomize_initial_coordinates = lround(params.parameters["randomize_coord"]);

  // Output filenames for different times
  ofstream fout_filenames(params.output_prefix + ".filenames" + ".dat");

  double tbeg = LennardJones::get_wall_time();

  // Loop over the events
  for(int n_event = 0; n_event < nevents; ++n_event) {

    double tbegevt = LennardJones::get_wall_time();

    //cout << "Event # " << n_event + 1 << endl;

    syst.Reinitialize(config);
    if (!config.canonical) {
      while (!syst.RenormalizeVelocitiesToEnergy(ust))
        syst.Reinitialize(config);
    }

    //cout << setw(tabsize) << "t" << " " << setw(tabsize) << "u*" << " " << setw(tabsize) << "T*" << endl;

    double t = 0.;
    // Output the initial state
    RunBoxEventsFunctions::OutputEvent(syst,
                                       params.output_prefix + ".t" + to_string(t) + ".dat",
                                       n_event + 1);
    double tout_next = t + dt_out;

//    cout << setw(tabsize) << t << " " << setw(tabsize) << syst.U / syst.m_config.N << " " << setw(tabsize) << syst.T
//         << endl;
//    cout.flush();

    // Simulation
    while (t < tfin) {
      if (tout_next - t <= dt) {
        syst.Integrate(tout_next - t);
        t = tout_next;
        RunBoxEventsFunctions::OutputEvent(syst,
                                           params.output_prefix + ".t" + to_string(t) + ".dat",
                                           n_event + 1);
        tout_next += dt_out;

        if (n_event == 0) {
          fout_filenames << setw(tabsize) << t << "   " << params.output_prefix + ".t" + to_string(t) + ".dat" << endl;
        }

//        cout << setw(tabsize) << t << " " << setw(tabsize) << syst.U / syst.m_config.N << " " << setw(tabsize) << syst.T
//             << endl;
//        cout.flush();
        cout << "\r" <<
        "Event # " << n_event + 1 << " : t = " << t << " / " << tfin << " = " << t / tfin * 100 << "%\t";
        cout.flush();
      } else {
        syst.Integrate(dt);
        t += dt;
      }
    }

    double tendevt = LennardJones::get_wall_time();
    cout << "\r" << "Event # " << n_event + 1 << " finished in " << tendevt - tbegevt << " seconds" << endl;
    cout.flush();

    if (n_event == 0) {
      fout_filenames.close();
    }
  }

  double tend = LennardJones::get_wall_time();
  // Output time per event
  cout << "Time per event: " << (tend - tbeg) / nevents << " seconds" << endl;

  return 0;
}
 