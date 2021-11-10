#include "include/run-fluctuations-aux.h"
#include "NumberStatistics.h"
#include "cuda_runtime.h"

using namespace std;


int main(int argc, char *argv[])
{
	RunFluctuationsParameters params;

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
	if (!config.canonical) {
		syst.RenormalizeVelocitiesToEnergy(ust);
	}


	// Create counters
	double dalpha = params.parameters["subvolume_spacing"];
	RunFluctuationsFunctions::RDF_average RDF(5.0, 0.03);
	RunFluctuationsFunctions::CoordFlucsAverage flucs_coord_X(0, dalpha);
	RunFluctuationsFunctions::CoordFlucsAverage flucs_coord_Y(1, dalpha);
	RunFluctuationsFunctions::CoordFlucsAverage flucs_coord_Z(2, dalpha);
	RunFluctuationsFunctions::CoordFlucsAverage flucs_coord_cube(3, dalpha);
	RunFluctuationsFunctions::MomentumFlucsAverage flucs_velocity_z(syst, 2, dalpha);


	double t = 0.;

	// Equilibration phase
	while (t < teq) {
		syst.Integrate(dt);
		t += dt;
	}

	// Prepare output
	cout << setw(15) << "t*"   << " ";
	cout << setw(15) << "u*"   << " ";
	cout << setw(15) << "T*"   << " ";
	cout << setw(15) << "Z"    << " ";
	cout << setw(15) << "<u*>" << " ";
	cout << setw(15) << "<T*>" << " ";
	cout << setw(15) << "<Z>"  << " ";
	cout << setw(15) << "<w>/(1-x)" << " ";
	cout << endl;

	ofstream fout(params.output_prefix + string(".TimeDep.txt"));
	fout << setw(15) << "t*" << " ";
	fout << setw(15) << "u*" << " ";
	fout << setw(15) << "T*" << " ";
	fout << setw(15) << "Z" << " ";
	fout << setw(15) << "<u*>" << " ";
	fout << setw(15) << "<T*>" << " ";
	fout << setw(15) << "<Z>" << " ";
	fout << setw(15) << "<wx>/(1-x)" << " ";
	fout << setw(15) << "<wy>/(1-x)" << " ";
	fout << setw(15) << "<wz>/(1-x)" << " ";
	fout << setw(15) << "<vx>" << " ";
	fout << setw(15) << "<vy>" << " ";
	fout << setw(15) << "<vz>" << " ";
	fout << endl;

	// Production phase
	double tot_N = 0.0, tot_N2 = 0.0;  // Fluctuations at \alpha = 0.5
	int totIters = 0;  // Total iterations
	syst.resetAveraging();
	while (t < tfin || tfin < 0.) {
		syst.Integrate(dt);
		t += dt;
		totIters++;

		// Update statistics
		{
			RDF.AddTimeStep(syst);
			flucs_coord_X.AddTimeStep(syst);
			flucs_coord_Y.AddTimeStep(syst);
			flucs_coord_Z.AddTimeStep(syst);
			flucs_coord_cube.AddTimeStep(syst);
			flucs_velocity_z.AddTimeStep(syst);
		}

		int tN = RunFluctuationsFunctions::GetNSubsystem(syst, 0.5, 2);
		tot_N += tN;
		tot_N2 += tN * tN;

		// Output
		if (totIters % 1000 == 0) {
			cout << setw(15) << t << " ";
			cout << setw(15) << syst.U / syst.m_config.N << " ";
			cout << setw(15) << syst.T << " ";
			cout << setw(15) << syst.P / (syst.m_config.rho * syst.T) << " ";
			cout << setw(15) << syst.av_U_tot / syst.av_iters / syst.m_config.N << " ";
			double Tav = syst.av_T_tot / syst.av_iters;
			cout << setw(15) << Tav << " ";
			cout << setw(15) << syst.av_p_tot / syst.av_iters / (syst.m_config.rho * Tav) << " ";

			double Nav = tot_N / totIters;
			double N2av = tot_N2 / totIters;


			cout << setw(15) << (N2av - Nav * Nav) / Nav / (1. - 0.5) << " ";
			cout << endl;


			fout << setw(15) << t << " ";
			fout << setw(15) << syst.U / syst.m_config.N << " ";
			fout << setw(15) << syst.T << " ";
			fout << setw(15) << syst.P / (syst.m_config.rho * syst.T) << " ";
			fout << setw(15) << syst.av_U_tot / syst.av_iters / syst.m_config.N << " ";
			fout << setw(15) << Tav << " ";
			fout << setw(15) << syst.av_p_tot / syst.av_iters / (syst.m_config.rho * Tav) << " ";
			//fout << setw(15) << (N2av - Nav * Nav) / Nav / (1. - 0.5) << " ";
			{
				RunFluctuationsFunctions::CoordFlucsAverage *flucs;
				for(int ic = 0; ic < 3; ++ic) {
					if (ic == 0)
						flucs = &flucs_coord_X;
					else if (ic == 1)
						flucs = &flucs_coord_Y;
					else
						flucs = &flucs_coord_Z;

					int ti = (*flucs).alphas.size() / 2;
					double alpha = (*flucs).alphas[ti];
					double tNav  = (*flucs).totsN[ti]  / (*flucs).iters;
					double tN2av = (*flucs).totsN2[ti] / (*flucs).iters;
					double tw = (tN2av - tNav * tNav) / tNav;
					fout << setw(15) << tw / (1. - alpha) << " ";
				}
			}


			auto velo = RunFluctuationsFunctions::GetAvVel(syst);
			fout << setw(15) << velo[0] << " ";
			fout << setw(15) << velo[1] << " ";
			fout << setw(15) << velo[2] << " ";
			fout << endl;
			fout.flush();

			// Other outputs
			RDF.PrintToFile(params.output_prefix + string(".RDF.dat"));
			flucs_coord_X.PrintToFile(params.output_prefix + string(".flucsX.dat"));
			flucs_coord_Y.PrintToFile(params.output_prefix + string(".flucsY.dat"));
			flucs_coord_Z.PrintToFile(params.output_prefix + string(".flucsZ.dat"));
			flucs_coord_cube.PrintToFile(params.output_prefix + string(".flucsCube.dat"));
			flucs_velocity_z.PrintToFile(params.output_prefix + string(".flucsVz.dat"));
		}
	}

  return 0;
}
 