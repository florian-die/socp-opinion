#include "socp/shooting.hpp"
#include "opinion_functions.hpp"
#include "opinion_model.hpp"

#include <vector>
#include <iostream>
#include <cmath>
#include <time.h>

// main function
int main(int argc, char** argv)
{
	std::cout << std::endl;
	std::cout << "----- SOCP for Opinion Dynamics -----" << std::endl;
	std::cout << std::endl;

	/* -------------- Parameters ---------------------------------------------- */
	int N = 6; // number agents (leader excluded)

	opinion my_opinion(N);

	my_opinion.model_params.N = N;
	my_opinion.model_params.sigma = 5.0;

	double eta = 0.5; // final neigborhood
	my_opinion.problem_params.eta = eta;

	my_opinion.homotopy_params.u = 1.0;
	my_opinion.homotopy_params.h = 0.0;

	my_opinion.ode_params.steps = 100;

	my_opinion.output_params.file_name = "opinion_data.dat";

	shooting my_shooting(my_opinion, 1, 1);

	/* -------------- Final Mode ---------------------------------------------- */
	// final time mode
	int mode_tf = opinion::MODE::FREE;

	// opinion final mode
	std::vector<int> mode_Yf(N+1,opinion::MODE::FREE);
	mode_Yf[1] = opinion::MODE::FIXED;
	mode_Yf[N] = opinion::MODE::FIXED;

	my_shooting.SetMode(mode_tf,mode_Yf);

	/* -------------- Intial & Final States ----------------------------------- */

	// initial & final times
	real ti = 0; // fixed
	real tf = 17.0; // free

	// initial opinions (fixed)
	opinion::mstate Yi(N+1);
	Yi[0] = 0.0;
	Yi[1] =	3.0;
	Yi[2] =	3.5;
	Yi[3] =	4.0;
	Yi[4] =	12.0;
	Yi[5] =	13.0;
	Yi[6] =	14.0;

	// final opinions
	opinion::mstate Yf(N+1,0.0);
	// Yf[0] = 5.4; // free
	Yf[1] =	-eta; // fixed
	// Yf[2] =	-0.1; // free
	// Yf[3] =	-0.1; // free
	// Yf[4] =	-0.1; // free
	// Yf[5] =	-0.1; // free
	Yf[N] =	eta; // fixed

	// initial costates (guess)
	opinion::mstate Pi(N+1,0.0);
	Pi[1] = -0.5;

	// double a = -0.5;
	// double b = - Pi[1] - opinion_functions::h(Yi[N]);
	// double c = - Pi[1]*Pi[1] - Pi[1]*opinion_functions::h(Yi[1])+1;
	// double d = b*b-4*a*c;
	//
	// Pi[N] = (-b-sqrt(d))/2/a;
	// std::cout << "Pi[N] = " << Pi[N] << std::endl;

	Pi[N] = 2.0;

	// final costates (not used)
	opinion::mstate Pf(N+1,0.0);

	// concatenation
	opinion::mstate Xi = my_opinion.Fuse(Yi,Pi);
	opinion::mstate Xf = my_opinion.Fuse(Yf,Pf);

	my_shooting.InitShooting(ti, Xi, tf, Xf);

	/* -------------- Solving initial problem---------------------------------- */

	std::cout << "1) Solving initial OCP... " << std::endl;
	double time1 = clock();
	int info = my_shooting.SolveOCP();
	double time2 = clock();
	double time = (double)(time2 - time1) / CLOCKS_PER_SEC;
	std::cout << "  - Algo returned " << info << std::endl;
	std::cout << "  - Computing time : " << time << " sec" << std::endl;

	if (info != 1)
	{
		std::cout << "!! Optimisation failed !!" << std::endl;
		return info;
	}

	/* -------------- Continuation on quadratic control cost ------------------ */

	std::cout << "2) Continuation on parameter u... " << std::endl;
	time1 = clock();
	info = my_shooting.SolveOCP(0.1, my_opinion.homotopy_params.u, 0.6);
	time2 = clock();
	time = (time2 - time1) / CLOCKS_PER_SEC;
	std::cout << "  - Algo returned " << info << std::endl;
	std::cout << "  - Computing time : " << time << " sec" << std::endl;

	if (info != 1)
	{
		std::cout << "!! Optimisation failed !!" << std::endl;
		return info;
	}

	/* -------------- Write solution in a file -------------------------------- */
	my_shooting.Trace();

	return 0;
}
