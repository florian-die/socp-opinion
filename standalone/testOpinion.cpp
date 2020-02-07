/*
 * testDoubleIntegrator.cpp
 *
 *  Created on: January 30, 2018
 *      Author: Bruno HERISSE
 */
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "socp/shooting.hpp"
#include "opinion.hpp"

// main function
int main(int argc, char** argv) {

	/********************************************************/
	/***************** OCP initialization  ******************/
	/********************************************************/

	std::cout << "Initialization... ";
	/// new map object (environment parameters, obstacles, etc.)
	//map* my_map = new map();
	/// new model object
	std::string my_fileTrace("opinion_data.dat");
	opinion* my_opinion = new opinion(NULL, my_fileTrace);

	// new shooting object
	shooting* my_shooting = new shooting(*my_opinion, 1, 1);	// model / number for multiple shooting (1 if simple shooting)
	int mode_tf = 1;											// 0 for fixed final time, 1 for free final time
	std::vector<int> mode_X;
	mode_X = std::vector<int>(my_opinion->GetDim(),0);
	for (int k=2;k<my_opinion->GetDim()-1; k++){
		mode_X[k] = 1;
	}
	mode_X[0] = 1; mode_X[1] = 0; mode_X[my_opinion->GetDim()-1] = 0;
	my_shooting->SetMode(mode_tf,mode_X);							// free final time / static final state (default)

	// set initial guess for the shooting
	real ti = 0;											// initial guess for ti
	model::mstate Xi(2*my_opinion->GetDim());		// GetDim return the dimension of the dynamic model. With the adjoint vector, the full state is twice this dimension

	//Xi[0] =		0.0;				//x
	//Xi[1] =		1;					//y
	//Xi[my_opinion->GetDim()-1] =		6;				//z
	//Xi[my_opinion->GetDim()] =		    0.01;				//p_x					// the adjoint vector is the guess
	//Xi[my_opinion->GetDim()+1] =		0.01;				//p_y
	//Xi[2*my_opinion->GetDim()-1] =		0.01;				//p_z
	//for (int k=2;k<my_opinion->GetDim()-1; k++){
	//	Xi[k] = 1 + 5.0*(k-1)/10;
	//	Xi[k + my_opinion->GetDim()] =		0.01;
	//}

	Xi[0] =		0.0;				//x
	Xi[1] =		3;					//y
	Xi[2] =		3.5;
	Xi[3] =		4;
	Xi[4] =		12;					//y
	Xi[5] =		13;
	Xi[6] =		14;

	//Xi[my_opinion->GetDim()-1] =		6;				//z
	Xi[my_opinion->GetDim()] =		    0.01;				//p_x					// the adjoint vector is the guess
	Xi[my_opinion->GetDim()+1] =		0.01;				//p_y
	Xi[2*my_opinion->GetDim()-1] =		0.01;				//p_z
	for (int k=2;k<my_opinion->GetDim()-1; k++){
		//Xi[k] = 1 + 5.0*(k-1)/10;
		Xi[k + my_opinion->GetDim()] =		0.01;
	}

	real tf = 5;											// initial guess for tf

	model::mstate Xf(2*my_opinion->GetDim());				// No need to specify the adjoint vector for Xf
	Xf[0] =		1.0;				//x
	Xf[1] =		-0.5;				//y
	Xf[my_opinion->GetDim()-1] =		0.5;				//z
	for (int k=2;k<my_opinion->GetDim()-1; k++){
		Xf[k] = 0.0;
	}

	my_shooting->InitShooting(ti, Xi, tf, Xf);					// init shooting with this guess
	std::cout << "done" << std::endl;

	/********************************************************/
	/*************** compute optimal control ****************/
	/********************************************************/
	// to return computing time
	long int time1, time2;
	double time;

	// if you want to use a different desired state from the intial guess
	//my_shooting->SetDesiredState(ti, Xi, tf, Xf);

	// print initial parameters for shooting (=tf+adjoint)
	std::cout << "	INIT : ";
	if (mode_tf==1){
		std::cout << "tf = " << my_shooting->GetParameters(mode_tf + (my_opinion->GetDim()-1)) << ", ";
	}else{
		std::cout << "tf = " << tf << ", ";
	}
	std::cout << "p_x = " << my_shooting->GetParameters(0) << ", ";
	std::cout << "p_y = " << my_shooting->GetParameters(1) << ", ";
	std::cout << "p_z = " << my_shooting->GetParameters(my_opinion->GetDim()-1) << std::endl;

	// solve OCP
	std::cout << "Solve OCP... ";
	time1 = clock();
	int info = my_shooting->SolveOCP(0.0);		// solve OCP with continuation step 0
	time2 = clock();
	time = (double)(time2 - time1) / CLOCKS_PER_SEC;
	std::cout << "Algo returned " << info << ", ";
	std::cout << "Computing time : " << time << " secondes." << std::endl;

	// print solution
	std::cout << "	SOL :  ";
	if (mode_tf==1){
		std::cout << "tf = " << my_shooting->GetParameters(mode_tf + (my_opinion->GetDim()-1)) << ", ";
	}else{
		std::cout << "tf = " << tf << ", ";
	}
	std::cout << "p_x = " << my_shooting->GetParameters(0) << ", ";
	std::cout << "p_y = " << my_shooting->GetParameters(1) << ", ";
	std::cout << "p_z = " << my_shooting->GetParameters(my_opinion->GetDim()-1) << std::endl;

	/********************************************************/
	/*********** continuation on a model parameter **********/
	/********************************************************/
	// model parameters can be changed using continuation
	std::cout << "Continuation on parameter muh... " << std::endl;
	if (info==1)	info = my_shooting->SolveOCP(0.1, my_opinion->GetParameterData()->muh, 0.2);

	// std::cout << "Continuation on parameter muU... " << std::endl;
	// //if (info==1)	info = my_shooting->SolveOCP(0.1, my_opinion->GetParameterData()->muU, 0.04);
	// if (info==1)	info = my_shooting->SolveOCP(0.1, my_opinion->GetParameterData()->muU, 0.04);
	// std::cout << "OK = " << info << ", done" << std::endl;

	// print solution
	std::cout << "	SOL :  ";
	if (mode_tf==1){
		std::cout << "tf = " << my_shooting->GetParameters(mode_tf + (my_opinion->GetDim()-1)) << ", ";
	}else{
		std::cout << "tf = " << tf << ", ";
	}
	std::cout << "p_x = " << my_shooting->GetParameters(0) << ", ";
	std::cout << "p_y = " << my_shooting->GetParameters(1) << ", ";
	std::cout << "p_z = " << my_shooting->GetParameters(my_opinion->GetDim()-1) << std::endl;

	/********************************************************/
	/*************** plot management ************************/
	/********************************************************/
	// print the trajectory in my_fileTrace from the last obtained solution
	my_shooting->Trace(); 	// trace
	// call plot scripts (this can be changed)
	#ifdef WIN32 /* Windows */
		//system("matlab -minimize -r \"run('../../trace/opinion/mplot.m');exit;\""); //-nodisplay -nosplash -nodesktop
	#elif defined (linux) /* Linux */
		// system("gnuplot gplot.sh");
		system("python trace_opinion_full.py");
	#else /* non supported plateform */
		#error not defined for this platform
	#endif

	//system("pause");

	return 0;
}
