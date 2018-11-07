// ##################################################
// #   main.cpp - finally revised on Jun 2018     	#
// #   coded by Jorge Escribano                     #
// #   Copyright (C) 2016-2018, Jorge Escribano,    #
// #   All rights reserved.                         #
// ##################################################

// Main file 

#include "class/Monolayer.hpp"



int main (int argc, char**argv)
{
	if (argc < 2) 
	{       
		fprintf(stderr, "Falta el indice.\n");   // Stop if error
		//exit(1);
		argv[1]=0;
	}
	char *index_cluster = argv[1];  // Gets hermes index


	Monolayer System;
	int ctFreq=0,ctParaview=0;
	System.index_cluster_Aux=index_cluster;






/*
  	std::chrono::duration<double> time_loop, time_loop2, time_loop3;
  	std::chrono::high_resolution_clock::time_point start, end, begin;
  	start = std::chrono::high_resolution_clock::now();
  	begin = std::chrono::high_resolution_clock::now();
	end = std::chrono::high_resolution_clock::now();
	time_loop= std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
	time_loop2= std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
	time_loop3= std::chrono::duration_cast<std::chrono::duration<double>>(end - start);

*/


	//System.BalanceToy();

	for (int i=0; i< System.GETTimeSlots(); ++i)
	{
		std::cout << "Step:  " << i << std::endl;

//		start = std::chrono::high_resolution_clock::now();
		System.Binding();
//		end = std::chrono::high_resolution_clock::now();
//		time_loop3 += std::chrono::duration_cast<std::chrono::duration<double>>(end - start);

		System.CheckReBinding();
		System.ForceGeneration();

		System.CellRemodel();
		
//		start = std::chrono::high_resolution_clock::now();
		System.BalanceLangevin();
//		end = std::chrono::high_resolution_clock::now();
//		time_loop += std::chrono::duration_cast<std::chrono::duration<double>>(end - start);

		System.UpdateStep();
		System.CheckUnbinding();
		System.UpdateCadhProp();	
		System.UpdateAdhPointCts();

//		start = std::chrono::high_resolution_clock::now();
		System.CalculateMembraneDistance();
//		end = std::chrono::high_resolution_clock::now();
//		time_loop2 += std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
//		System.printVTK();


		if ( System.GETprintFilesVTK() )
		{	 
			if (  ctFreq >= System.GETplotFreq() )
			{
				ctParaview++;
				System.printVTK(ctParaview);
				
				ctFreq=0;
			}
			ctFreq++;
		}
		
		
		System.UpdateGaps(i);
		std::cout << "Update: OK:  " << std::endl;
		System.CheckGapFormation(i);
		std::cout << "Formation OK:  " << std::endl;
		System.PrintGAPStep(index_cluster);
		System.PrintStepStat(index_cluster);

	//	System.printDataStepVertex(index_cluster,i);
	}



	//System.printDataGlobalVertex(index_cluster);
	//System.PrintGAP(index_cluster);





	System.PrintGlobalStat(index_cluster);


/*
	end = std::chrono::high_resolution_clock::now();

	std::cout << "time substrate_balance: " << time_loop.count()<< std::endl;
	std::cout << "percentage Langevin: " << time_loop.count()/std::chrono::duration_cast<std::chrono::duration<double>>(end-begin).count()<< std::endl;
	std::cout << "percentage UdateGap: " << time_loop2.count()/std::chrono::duration_cast<std::chrono::duration<double>>(end-begin).count()<< std::endl;
	std::cout << "percentage Binding: " << time_loop3.count()/std::chrono::duration_cast<std::chrono::duration<double>>(end-begin).count()<< std::endl;
	//std::cout << "loop" << duration_cast<milliseconds>(stop - start).count()<< std::endl;
         
*/

	System.matlabPlot(index_cluster);

//	system ("rm -r results");
	system ("mkdir results");
    system ("mv *.vtk ./results/");

    system ("mv *.txt ./results/");
    system ("cp *.eps ./results/");
    system ("cp *.conf ./results/");
    system ("cp *.jpg ./results/");
    system ("cp *.m ./results/");


	return 0;
}