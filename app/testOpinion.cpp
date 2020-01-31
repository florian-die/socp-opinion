#include "socp/shooting.hpp"
#include "opinion_model.hpp"

#include <vector>
#include <iostream>

// main function
int main(int argc, char** argv)
{
	/* -------------- Parameters ---------------------------------------------- */
	int N = 6; // number agents (leader excluded)

	opinion my_opinion(N);

	my_opinion.model_params.N = N;
	my_opinion.model_params.sigma = 5.0;

	double eta = 0.5; // final neigborhood
	my_opinion.problem_params.eta = eta;

	my_opinion.homotopy_params.u = 1.0;
	my_opinion.homotopy_params.h = 0.0;

	my_opinion.ode_params.steps = 500;

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
	real tf = 20.0; // free

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
	Yf[0] = 12.0; // free
	Yf[1] =	eta; // fixed
	Yf[N] =	-eta; // fixed

	// initial costates (guess)
	opinion::mstate Pi(N+1,0.0);
	Pi[1] = 0.02;
	Pi[N] = -0.01;

	// final costates (not used)
	opinion::mstate Pf(N+1,0.0);

	// concatenation
	opinion::mstate Xi = my_opinion.Fuse(Yi,Pi);
	opinion::mstate Xf = my_opinion.Fuse(Yf,Pf);

	my_shooting.InitShooting(ti, Xi, tf, Xf);

	/* -------------- Solving initial problem---------------------------------- */

	std::cout << "Solve OCP... ";
	long int time1 = clock();
	int info = my_shooting.SolveOCP(0.0);
	long int time2 = clock();
	double time = (double)(time2 - time1) / CLOCKS_PER_SEC;
	std::cout << "Algo returned " << info << ", ";
	std::cout << "Computing time : " << time << " secondes." << std::endl;


	/* -------------- Write solution in a file -------------------------------- */
	my_shooting.Trace();

	return 0;
}
