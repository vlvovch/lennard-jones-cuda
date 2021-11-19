#include "include/run-isotherm-aux.h"
#include "../run-fluctuations/include/run-fluctuations-aux.h"
#include "NumberStatistics.h"
#include "cuda_runtime.h"

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
		cudaGetDeviceCount(&deviceCount);

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
	cout << setw(15) << "<p*>" << " ";
	cout << setw(15) << "d<p*>" << " ";
	cout << setw(15) << "<Z>" << " ";
	cout << setw(15) << "d<Z>" << " ";
	cout << endl;


	ofstream fout(params.output_prefix + string(".dat"));
	fout << setw(15) << "rho*" << " ";
	fout << setw(15) << "t*" << " ";
	fout << setw(15) << "T*" << " ";
	fout << setw(15) << "<u*>" << " ";
	fout << setw(15) << "d<u*>" << " ";
	fout << setw(15) << "<p*>" << " ";
	fout << setw(15) << "d<p*>" << " ";
	fout << setw(15) << "<Z>" << " ";
	fout << setw(15) << "d<Z>" << " ";
	fout << setw(15) << "<h>" << " ";
	fout << setw(15) << "d<h>" << " ";
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

		TimeAverage uav, Pav, Nav, Pxyav;
		int iter_dec_max = 200;
		int diter_dec = 1;
		AnalyzeExpDecay Pav_dec(diter_dec, iter_dec_max);
		AnalyzeExpDecay Pav_dec_alt(diter_dec, iter_dec_max);
		AnalyzeExpDecay Pxyav_dec(diter_dec, iter_dec_max);
		AnalyzeExpDecay Pxyav_dec_alt(diter_dec, iter_dec_max);
		AnalyzeExpDecay Nav_dec(diter_dec, iter_dec_max*10);
		std::vector<double> Ps;
		std::vector<double> Pxys;
		std::vector<double> Ps_roll(iter_dec_max + 1);
		std::vector<double> Pxys_roll(iter_dec_max + 1);
		std::vector<double> Ns;

		SampleMoments::NumberStatistics uloc, Ploc, utot, Ptot;
		SampleMoments::NumberStatistics hloc,  htot; // Enthalpy

		double t_next_eq = t + teq;

		while (t < tfin) {
			uav.AddObservation(syst.U / syst.m_config.N);
			Pav.AddObservation(syst.P);
			Pxyav.AddObservation(syst.Pshear);

			Ps.push_back(syst.P);
			Pav_dec.AddObservation(Ps);

			Pxys.push_back(syst.Pshear);
			Pxyav_dec.AddObservation(Pxys);
			double T = syst.T;

			int tN = RunFluctuationsFunctions::GetNSubsystem(syst, 0.5, 2);
			Ns.push_back(tN);
			Nav.AddObservation(tN);
			Nav_dec.AddObservation(Ns);

			if (totIters != 0 && totIters % iter_dec_max == 0) {
				Ps_roll[iter_dec_max] = syst.P;
				Pxys_roll[iter_dec_max] = syst.Pshear;
				Pav_dec_alt.AddSegment(Ps_roll);
				Pxyav_dec_alt.AddSegment(Pxys_roll);
			}
			Ps_roll[totIters % iter_dec_max] = syst.P;
			Ps_roll[totIters % iter_dec_max] = syst.Pshear;

			uloc.AddObservation(syst.U / syst.m_config.N);
			Ploc.AddObservation(syst.P);
			hloc.AddObservation(syst.U * syst.m_config.rho / syst.m_config.N + syst.P);

			if (t > t_next_eq) {
				utot.AddObservation(uloc.GetMean());
				Ptot.AddObservation(Ploc.GetMean());
				htot.AddObservation(hloc.GetMean());

				uloc.Clear();
				Ploc.Clear();
				hloc.Clear();

				t_next_eq += teq;
			}

			totIters++;
			if (totIters % 1000 == 0) {
				cout << setw(15) << rho << " ";
				cout << setw(15) << t << " ";
				cout << setw(15) << T << " ";
				cout << setw(15) << utot.GetMean() << " ";
				cout << setw(15) << utot.GetMeanError() << " ";
				cout << setw(15) << Ptot.GetMean() << " ";
				cout << setw(15) << Ptot.GetMeanError() << " ";
				cout << setw(15) << Ptot.GetMean() / rho / syst.m_config.T0 << " ";
				cout << setw(15) << Ptot.GetMeanError() / rho / syst.m_config.T0 << " ";

				cout << endl;

				Pav_dec.PrintToFile(params.output_prefix + ".PressureDecay.rho." + std::to_string(rho) + ".dat", 
					dt, 
					0.5 * Pav.GetS() * dt,
					syst.m_config.N / syst.m_config.rho / syst.m_config.T0
				);
				Pav_dec_alt.PrintToFile(params.output_prefix + ".PressureDecayIndep.rho." + std::to_string(rho) + ".dat", 
					dt, 
					0.5 * Pav.GetS() * dt,
					syst.m_config.N / syst.m_config.rho / syst.m_config.T0
				);
				Nav_dec.PrintToFile(params.output_prefix + ".NumberDecay.rho." + std::to_string(rho) + ".dat", 
					dt, 
					0.5 * Nav.GetS() * dt,
					2. / Nav_dec.mean.GetVariance()
					);
				Pxyav_dec.PrintToFile(params.output_prefix + ".PxyDecay.rho." + std::to_string(rho) + ".dat", 
					dt, 
					0.5 * Pxyav.GetS() * dt,
					syst.m_config.N / syst.m_config.rho / syst.m_config.T0
				);
				Pxyav_dec_alt.PrintToFile(params.output_prefix + ".PxyDecayIndep.rho." + std::to_string(rho) + ".dat",
					dt,
					0.5 * Pav.GetS() * dt,
					syst.m_config.N / syst.m_config.rho / syst.m_config.T0
				);
			}

			syst.Integrate(dt);
			t += dt;
		}

		fout << setw(15) << rho << " ";
		fout << setw(15) << t - dt << " ";
		fout << setw(15) << syst.T << " ";
		fout << setw(15) << utot.GetMean() << " ";
		fout << setw(15) << utot.GetMeanError() << " ";
		fout << setw(15) << Ptot.GetMean() << " ";
		fout << setw(15) << Ptot.GetMeanError() << " ";
		fout << setw(15) << Ptot.GetMean() / rho / syst.m_config.T0 << " ";
		fout << setw(15) << Ptot.GetMeanError() / rho / syst.m_config.T0 << " ";
		fout << setw(15) << htot.GetMean() << " ";
		fout << setw(15) << htot.GetMeanError() << " ";

		double Pcur = Ptot.GetMean();
		double Pcurerr = Ptot.GetMeanError();

		double tdrho = drho;
		if (rho == rhomin)
			tdrho = rhomin;

		double w = syst.m_config.T0 / ((Pcur - Pprev) / tdrho);
		double werr = syst.m_config.T0 * tdrho * sqrt(Pcurerr * Pcurerr + Ppreverr * Ppreverr) / (Pcur - Pprev) / (Pcur - Pprev);

		fout << setw(15) << w << " ";
		fout << setw(15) << werr << " ";

		Pprev = Pcur;
		Ppreverr = Pcurerr;

		fout << endl;

		fout.flush();

		Pav_dec.PrintToFile(params.output_prefix + ".PressureDecay.rho." + std::to_string(rho) + ".dat",
			dt,
			0.5 * Pav.GetS() * dt,
			syst.m_config.N / syst.m_config.rho / syst.m_config.T0
		);
		Pav_dec_alt.PrintToFile(params.output_prefix + ".PressureDecayIndep.rho." + std::to_string(rho) + ".dat",
			dt,
			0.5 * Pav.GetS() * dt,
			syst.m_config.N / syst.m_config.rho / syst.m_config.T0
		);
		Nav_dec.PrintToFile(params.output_prefix + ".NumberDecay.rho." + std::to_string(rho) + ".dat",
			dt,
			0.5 * Nav.GetS() * dt,
			2. / Nav_dec.mean.GetVariance()
		);
		Pxyav_dec.PrintToFile(params.output_prefix + ".PxyDecay.rho." + std::to_string(rho) + ".dat",
			dt,
			0.5 * Pxyav.GetS() * dt,
			syst.m_config.N / syst.m_config.rho / syst.m_config.T0
		);
		Pxyav_dec_alt.PrintToFile(params.output_prefix + ".PxyDecayIndep.rho." + std::to_string(rho) + ".dat",
			dt,
			0.5 * Pav.GetS() * dt,
			syst.m_config.N / syst.m_config.rho / syst.m_config.T0
		);
	}


	return 0;
}
