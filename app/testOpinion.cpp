#include "socp/shooting.hpp"
#include "opinion_model.hpp"

#include <vector>

// main function
int main(int argc, char** argv)
{
	int N = 6;

	opinion my_opinion(N);

	my_opinion.model_params.N = N;
	my_opinion.model_params.sigma = 5.0;

	double eta = 0.5;
	my_opinion.problem_params.eta = eta;

	my_opinion.homotopy_params.u = 1.0;
	my_opinion.homotopy_params.h = 0.0;

	my_opinion.ode_params.steps = 100;

	shooting my_shooting(my_opinion, 1, 1);

	int mode_tf = 1; // 0 fixed / 1 free

	std::vector<int> mode_Xf = std::vector<int>(2*N+2,0);
	mode_Xf[0] = 1;

	for (int i = 2; i < N-1; i++)
	{
		mode_Xf[i] = 1;
	}

	my_shooting.SetMode(mode_tf,mode_Xf);

	real ti = 0;
	real tf = 5.0;

	opinion::mstate Yi(N+1);
	Yi[0] = 0;
	Yi[1] =	3.0;
	Yi[2] =	3.5;
	Yi[3] =	4.0;
	Yi[4] =	12.0;
	Yi[5] =	13.0;
	Yi[6] =	14.0;

	opinion::mstate Yf(N+1);
	Yf[0] = 0.0;
	Yf[1] =	eta;
	Yf[N] =	-eta;
	for (int i = 2; i < N; i++)
	{
		Yf[i] = 0.0;
	}

	opinion::mstate Pi(N+1);
	Pi[0] = 0.0;
	Pi[0] = 0.1;
	Pi[0] = 0.0;
	Pi[0] = 0.0;
	Pi[0] = 0.0;
	Pi[0] = 0.0;
	Pi[0] = 0.1;

	opinion::mstate Pf(N+1);
	for (int i = 0; i <= N; i++)
	{
		Pf[i] = 0.0;
	}

	opinion::mstate Xi = my_opinion.Fuse(Yi,Pi);
	opinion::mstate Xf = my_opinion.Fuse(Yf,Pf);

	my_shooting.InitShooting(ti, Xi, tf, Xf);

	my_shooting.SolveOCP(0.0);

	return 0;
}
