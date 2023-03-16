#include "include/run-isotherm-aux.h"
#include "NumberStatistics.h"
#ifdef USE_CUDA_TOOLKIT
#include "cuda_runtime.h"
#endif

using namespace std;


int main(int argc, char* argv[])
{
	RunIsothermParameters params;

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


	int N = lround(params.parameters["N"]);
	double Tst = params.parameters["T*"];

	// Density parameters
	double rhomin = params.parameters["rho*_min"];
	double rhomax = params.parameters["rho*_max"];
	double drho = params.parameters["drho*"];


	// Time parameters
	double dt = params.parameters["dt*"];
	double teq = params.parameters["teq"];
	double tfin = params.parameters["tfin"];

	// Prepare output
	cout << setw(15) << "rho*" << " ";
	cout << setw(15) << "t*" << " ";
	cout << setw(15) << "T*" << " ";
	cout << setw(15) << "<u*>" << " ";
	cout << setw(15) << "d<u*>" << " ";
	cout << setw(15) << "teq_u*" << " ";
	cout << setw(15) << "<p*>" << " ";
	cout << setw(15) << "d<p*>" << " ";
	cout << setw(15) << "teq_p*" << " ";
	cout << setw(15) << "<Z>" << " ";
	cout << setw(15) << "d<Z>" << " ";
	cout << endl;


	ofstream fout(params.output_prefix + string(".dat"));
	fout << setw(15) << "rho*" << " ";
	fout << setw(15) << "t*" << " ";
	fout << setw(15) << "T*" << " ";
	fout << setw(15) << "<u*>" << " ";
	fout << setw(15) << "d<u*>" << " ";
	fout << setw(15) << "teq_u*" << " ";
	fout << setw(15) << "<p*>" << " ";
	fout << setw(15) << "d<p*>" << " ";
	fout << setw(15) << "teq_p*" << " ";
	fout << setw(15) << "<Z>" << " ";
	fout << setw(15) << "d<Z>" << " ";
	fout << setw(15) << "w[N]" << " ";
	fout << setw(15) << "dw[N]" << " ";
	fout << endl;

	double Pprev = 0., Ppreverr = 0.;

	for (double rho = rhomin; rho <= rhomax; rho += drho) {
		cout << endl;

		MDSystem::MDSystemConfiguration config;
		config.N = N;
		config.T0 = Tst;
		config.rho = rho;
		config.useCUDA = lround(params.parameters["useCUDA"]);
		config.canonical = 1;

		MDSystem syst(config);
		syst.Reinitialize(config);

		double t = 0.;

		// Equilibration phase
		while (t < teq) {
			syst.Integrate(dt);
			t += dt;
		}

		int totIters = 0;

		TimeAverage uav, Pav;

		while (t < tfin) {
			uav.AddObservation(syst.U / syst.m_config.N);
			Pav.AddObservation(syst.P);
			double T = syst.T;


			totIters++;
			if (totIters % 1000 == 0) {
				cout << setw(15) << rho << " ";
				cout << setw(15) << t << " ";
				cout << setw(15) << T << " ";
				cout << setw(15) << uav.GetMean() << " ";
				cout << setw(15) << uav.GetMeanError() << " ";
				cout << setw(15) << uav.GetS() * dt << " ";
				cout << setw(15) << Pav.GetMean() << " ";
				cout << setw(15) << Pav.GetMeanError() << " ";
				cout << setw(15) << Pav.GetS() * dt << " ";
				cout << setw(15) << Pav.GetMean() / rho / syst.m_config.T0 << " ";
				cout << setw(15) << Pav.GetMeanError() / rho / syst.m_config.T0 << " ";
				cout << endl;
			}

			syst.Integrate(dt);
			t += dt;
		}

		fout << setw(15) << rho << " ";
		fout << setw(15) << t - dt << " ";
		fout << setw(15) << syst.T << " ";
		fout << setw(15) << uav.GetMean() << " ";
		fout << setw(15) << uav.GetMeanError() << " ";
		fout << setw(15) << uav.GetS() * dt << " ";
		fout << setw(15) << Pav.GetMean() << " ";
		fout << setw(15) << Pav.GetMeanError() << " ";
		fout << setw(15) << Pav.GetS() * dt << " ";
		fout << setw(15) << Pav.GetMean() / rho / syst.m_config.T0 << " ";
		fout << setw(15) << Pav.GetMeanError() / rho / syst.m_config.T0 << " ";

		double Pcur = Pav.GetMean();
		double Pcurerr = Pav.GetMeanError();

		double tdrho = drho;
		if (rho == rhomin)
			tdrho = rhomin;

		double w = syst.m_config.T0 / ((Pcur - Pprev) / tdrho);
		double werr = syst.m_config.T0 * tdrho * sqrt(Pcurerr * Pcurerr + Ppreverr * Ppreverr) / std::abs(Pcur - Pprev);

		fout << setw(15) << w << " ";
		fout << setw(15) << werr << " ";

		Pprev = Pcur;
		Ppreverr = Pcurerr;

		fout << endl;

		fout.flush();
	}


	return 0;
}
