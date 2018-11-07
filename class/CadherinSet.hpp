// ##################################################
// #   CadherinSet.hpp - finally revised on Jun 2018     #
// #   coded by Jorge Escribano                     #
// #   Copyright (C) 2016-2018, Jorge Escribano,    #
// #   All rights reserved.                         #
// ##################################################

#ifndef Cadherin_HPP
#define Cadherin_HPP

#include "../header.hpp"
#include "../config/SingleConfig.hpp"



struct Cadherin
{
	Cadherin(int first, int second, int lengthCtr, int cadherinNumberInit);
	Cadherin(int lengthCtr, int cadherinNumberInit);
	Cadherin(int lengthCtr);


	void UpdateCadhActivity();

	double force;
	double ubProb;
	double length;
	int firstAdhPoint;
	int secondAdhPoint;
	bool active;
	bool full;						// If cadherin 
	int cadherinNumber;
	
	int nGap;						// Number of times that a gap is generated
	bool gapCtAux;					// Only true at the step that it completely unbinds
	double gapTotTime;				// Total time withou adhesion
	std::vector<double> gapGenTime; // For each nGap

};

class CadherinSet
{
	public:
	 	 CadherinSet();
	 	 void UpdateCadhActivity();
	 	 void Binds(int ListID);
	 	 void UpdateUbProb(int ListID);
	 	 void Unbinds(int ListID);
 		 void updateCts(int timeSlots);
 		 void BindCadherin(int pointA, int pointB, int nInitAdh);
 		 void BindCadherinNew(int pointA, int pointB, int ListID, double distance);
 		 void AddFreeCadherin();



	 	// Getters
	 	 double GETbindLimit(){return bindLimit;};
	 	 double GETdensity(){return density;};
	 	 double GETEA(){return EA;};
	 	 double GETvisc(){return visc;};
	 	 double GETlengthBal(){return lengthBal;};
	 	 double GETKonRate(){return KonRate;};
	 	 double GETKBending(){return KBending;};
	 	 int GETclusterMaxNumber(){return clusterMaxNumber;};

	 	// Vector variables
	 	 std::vector<Cadherin> List;
    	 std::vector<int> ActiveCadh;	// Active is when it used for matrix
    	 std::vector<int> InactiveCadh;	// Inactive is when is not used
    	 std::vector<int> AdhCt;		// Adhesion cter	

	private:
		double density;
		double EA;
		double visc;
		double lengthBal;
		double bindLimit;
		double KonRate;
		double TimeStep;
		double KBending;
		// Catch bond variables
		double KCatch;
		double phiS;
		double phiC;
		double fAp;


		int clusterMaxNumber;


	
}; 




#endif /* Cadherin_HPP */