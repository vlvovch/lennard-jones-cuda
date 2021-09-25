#include "include/run-fluctuations-aux.h"
#include "NumberStatistics.h"
#include "cuda_runtime.h"

using namespace std;


int main(int argc, char *argv[])
{
	RunFluctuationsParameters params;

	// check CUDA
	{
		int deviceCount = 0;
		cudaGetDeviceCount(&deviceCount);

		if (deviceCount == 0 && lround(params.parameters["useCUDA"])) {
			std::cout << "Could not find a CUDA device! Calculations will be performed on CPU\n";
			params.parameters["useCUDA"] = 0;
		}
	}

	// Print all the run parameters
	params.OutputParameters();

	//cout << params.GetFullPrefix() << endl;

	int N = lround(params.parameters["N"]);
	double rhost = params.parameters["rho*"];
	double Tst = params.parameters["T*"];
	double ust = params.parameters["u*"];

  vector<double> fracs;
	{
		double dfr = params.parameters["subvolume_spacing"];
		for (double tfr = dfr; tfr <= 1.; tfr += dfr)
			fracs.push_back(tfr);
	}

	// Time parameters
	double dt = params.parameters["dt*"];
	double teq  = params.parameters["teq"];
	double tfin = params.parameters["tfin"];

	MDSystem::MDSystemConfiguration config;
	config.N = N;
	config.T0 = Tst;
	config.rho = rhost; 
	config.useCUDA = lround(params.parameters["useCUDA"]);
	config.canonical = lround(params.parameters["canonical"]);

	MDSystem syst(config);
	syst.Reinitialize(config);

	double t = 0.;


	// Equilibration phase
	while (t < teq) {
		syst.Integrate(dt);
		t += dt;
	}


	// Production phase
	double tot_N = 0.0, tot_N2 = 0.0;  // Fluctuations at \alpha = 0.5
	int totIters = 0;  // Total iterations
	syst.resetAveraging();
	while (t < tfin) {
		syst.Integrate(dt);
		t += dt;
		totIters++;

		int tN = RunFluctuationsParameters::GetNSubsystem(syst, 0.5, 2);
		tot_N += tN;
		tot_N2 += tN * tN;

		if (totIters % 100 == 0) {
			cout << setw(15) << t << " ";
			cout << setw(15) << syst.U / syst.m_config.N << " ";
			cout << setw(15) << syst.T << " ";
			cout << setw(15) << syst.P / (syst.m_config.rho * syst.m_config.T0) << " ";
			cout << setw(15) << syst.av_U_tot / syst.av_iters / syst.m_config.N << " ";
			cout << setw(15) << syst.av_T_tot / syst.av_iters << " ";
			cout << setw(15) << syst.av_p_tot / syst.av_iters / (syst.m_config.rho * syst.m_config.T0) << " ";

			double Nav = tot_N / totIters;
			double N2av = tot_N2 / totIters;


			cout << setw(15) << (N2av - Nav * Nav) / Nav / (1. - 0.5) << " ";
			cout << endl;
		}
	}

  return 0;
}
 