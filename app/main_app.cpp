#include "socp/shooting.hpp"
#include "opinion_functions.hpp"
#include "opinion_model.hpp"

#include <vector>
#include <iostream>
#include <cmath>
#include <time.h>
#include <stdlib.h>

// main function
int main(int argc, char** argv)
{
	std::cout << std::endl;
	std::cout << "----- SOCP for Opinion Dynamics -----" << std::endl;
	std::cout << std::endl;

	/* -------------- Parameters ---------------------------------------------- */
	int N = 6; // number agents (leader excluded)
	int M = N + 1;

	opinion my_opinion(N);

	my_opinion.model_params.N = N;
	my_opinion.model_params.sigma = 5.0;

	double eta = 0.5; // final neigborhood
	my_opinion.problem_params.eta = eta;

	my_opinion.solution_params.switching_times = std::vector<double>{ 0.0 };

	my_opinion.homotopy_params.u = 1.0;
	my_opinion.homotopy_params.h = 0.0;

	my_opinion.ode_params.steps = 1000;

	my_opinion.output_params.file_name = "opinion_data.dat";

	/* -------------- Initialisation ------------------------------------------ */

	int nMulti = 2;
	int nThread = 1;

	shooting my_shooting(my_opinion, nMulti, nThread);

	/* -------------- Final Mode ---------------------------------------------- */

	// final time mode
	int mode_tf = opinion::MODE::FREE;

	// opinion final mode
	std::vector<int> mode_Yf(M,opinion::MODE::FREE);
	mode_Yf[1] = opinion::MODE::FIXED;
	mode_Yf[N] = opinion::MODE::FIXED;

	my_shooting.SetMode(mode_tf,mode_Yf);

	/* -------------- Intial & Final States ----------------------------------- */

	// initial & final times
	real ti = 0; // fixed
	real tf = 17.0; // free

	// initial opinions (fixed)
	opinion::mstate Yi(M);
	Yi[0] = 0.0;
	Yi[1] =	3.0;
	Yi[2] =	3.5;
	Yi[3] =	4.0;
	Yi[4] =	12.0;
	Yi[5] =	13.0;
	Yi[6] =	14.0;

	// final opinions
	opinion::mstate Yf(M,0.0);
	Yf[1] =	-eta; // fixed
	Yf[N] =	eta; // fixed
	// others are free

	// initial costates (guess)
	opinion::mstate Pi(M,0.01);
	// Pi[1] = -0.5;
	// Pi[2] = -0.1;
	// Pi[3] = -0.1;
	// Pi[4] = 0.1;
	// Pi[5] = 0.1;
	// Pi[N] = 2.0;

	// final costates (not used)
	opinion::mstate Pf(M,0.0);

	// concatenation
	opinion::mstate Xi = my_opinion.Fuse(Yi,Pi);
	opinion::mstate Xf = my_opinion.Fuse(Yf,Pf);

	my_shooting.InitShooting(ti, Xi, tf, Xf);

	/* -------------- Init utilities ------------------------------------------ */

	double time1, time2, time;
	int info;

	/* -------------- Solving initial problem---------------------------------- */

	std::cout << "1) Solving initial OCP... " << std::endl;
	time1 = clock();
	info = my_shooting.SolveOCP();
	time2 = clock();
	time = (double)(time2 - time1) / CLOCKS_PER_SEC;
	std::cout << "  - Algo returned " << info << std::endl;
	std::cout << "  - Computing time : " << time << " sec"
						<< std::endl << std::endl;

	if (info != 1)
	{
		std::cout << "!! Optimisation failed !!" << std::endl;
		return info;
	}

	/* -------------- Continuation on interactions ----------------------------- */

	// std::cout << "2) Continuation on interactions... " << std::endl;
	// time1 = clock();
	// info = my_shooting.SolveOCP(0.01, my_opinion.homotopy_params.h, 0.3);
	// time2 = clock();
	// time = (time2 - time1) / CLOCKS_PER_SEC;
	// std::cout << "  - Algo returned " << info << std::endl;
	// std::cout << "  - Computing time : " << time << " sec"
	// 					<< std::endl << std::endl ;
	//
	// if (info != 1)
	// {
	// 	std::cout << "!! Optimisation failed !!" << std::endl;
	// 	return info;
	// }

	/* -------------- Continuation on quadratic control cost ------------------ */

	std::cout << "2) Continuation on quadratic control cost... " << std::endl;

	double step = 0.005;
	double end = 0.30;

	if (argc == 3)
	{
		step = strtod(argv[1],NULL);
		end = strtod(argv[2],NULL);
	}

	std::cout << "  - Continuation from " << my_opinion.homotopy_params.u
						<< " to " << end << " by " << step*100.0 << "% steps" << std::endl;

	time1 = clock();
	info = my_shooting.SolveOCP(step, my_opinion.homotopy_params.u, end);
	time2 = clock();
	time = (time2 - time1) / CLOCKS_PER_SEC;
	std::cout << "  - Algo returned " << info << std::endl;
	std::cout << "  - Computing time : " << time << " sec" << std::endl;
	std::cout << std::endl;

	if (info != 1)
	{
		std::cout << "!! Optimisation failed !!" << std::endl;
		return info;
	}

	/* -------------- Introducing switching times ------------------ */

	// // removing cost on control norm
	// my_opinion.homotopy_params.u = 0.0;
	//
	// // recover previous solution
	// std::vector<real> times(nMulti + 1);
	// std::vector<model::mstate> states(nMulti + 1);
	// my_shooting.GetSolution(times, states);
	// tf = times[1];
	//
	// // switching time (guessed from observation)
	// double t1 = 0.8;
	// my_opinion.solution_params.switching_times[0] = t1;
	// times[1] = t1;
	// times[2] = tf;
	//
	// states[0] = my_shooting.Move(times[0]);
	// states[1] = my_shooting.Move(times[1]);
	// states[2] = my_shooting.Move(times[2]);
	//
	// my_shooting.InitShooting(times, states);
	// // NB : set modes after init, not before
	//
	// // set time modes
	// std::vector<int> mode_t(nMulti+1, opinion::MODE::FREE);
	// mode_t[0] = opinion::MODE::FIXED;
	//
	// // set opinion modes
	// std::vector< std::vector<int> > mode_Y(nMulti+1);
	// mode_Y[0] = std::vector<int>(M, opinion::MODE::FIXED);
	// mode_Y[1] = std::vector<int>(M, 2);
	// mode_Y[nMulti] = mode_Yf;
	//
	// my_shooting.SetMode(mode_t, mode_Y);
	//
	// // solving
	// std::cout << "3) Introducing switching time" << std::endl;
	// time1 = clock();
	// info = my_shooting.SolveOCP();
	// time2 = clock();
	// time = (time2 - time1) / CLOCKS_PER_SEC;
	// std::cout << "  - Algo returned " << info << std::endl;
	// std::cout << "  - Computing time : " << time << " sec" << std::endl;
	//
	// if (info != 1)
	// {
	// 	std::cout << "!! Optimisation failed !!" << std::endl;
	// 	return info;
	// }
	//
	// std::cout << "  - Switching time : "
	// 					<< my_opinion.solution_params.switching_times[0]
	// 					<< " sec" << std::endl;
	//
	// my_shooting.GetSolution(times, states);
	//
	// std::cout << "  - Final time : " << times[2] << " sec" << std::endl;
	//
	// std::cout << "  - Final opinions : ["
	// 					<< states[2][1] << "," << states[2][N] << "]"
	// 					<< std::endl;

	/* -------------- Write solution in a file -------------------------------- */
	my_opinion.output_params.file_stream.open(my_opinion.output_params.file_name.c_str(), std::ios::trunc);
	my_opinion.output_params.file_stream.close();
	my_shooting.Trace();
	std::cout << "Written file : " << my_opinion.output_params.file_name.c_str() << std::endl;

	std::cout << std::endl;

	return 0;
}
