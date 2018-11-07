// ##################################################
// #   CadherinSet.cpp - finally revised on Jun 2018     	#
// #   coded by Jorge Escribano                     #
// #   Copyright (C) 2016-2018, Jorge Escribano,    #
// #   All rights reserved.                         #
// ##################################################

// Contain de adhesion complexes or cadherins class 

#include "CadherinSet.hpp"

CadherinSet::CadherinSet()
{
	int initialCadherinNumber, nCellRow, nInitAdh;
	double aristLength;
	int nAdhPointCell, i;

	// Read variables
	 SingleConfig DataConfig=SingleConfig::getConfig();
	 lengthBal=DataConfig.cadhLengthBal;
	 density=DataConfig.cadhDensity;
	 EA=DataConfig.cadhEA;
	 visc=DataConfig.cadhVisc;
	 bindLimit=DataConfig.bindLimitDistance;
	 KonRate=DataConfig.KonRate;
	 TimeStep=DataConfig.TimeStep;
	 clusterMaxNumber=DataConfig.clusterMaxNumber;
	 initialCadherinNumber=DataConfig.initialCadherinNumber;
	 nCellRow=DataConfig.nCell;
	 aristLength=DataConfig.aristLength;
	 KBending=DataConfig.KBendingCadh;

	 KCatch=DataConfig.KCatch;
	 phiS=DataConfig.phiS;
	 phiC=DataConfig.phiC;
	 fAp=DataConfig.fAp;


	 nAdhPointCell=floor((aristLength+0.1)/DataConfig.sepAdhPoint); 
	 aristLength=double(nAdhPointCell*DataConfig.sepAdhPoint);
	 nAdhPointCell=nAdhPointCell*6;


	 nInitAdh=(nCellRow-1)*5;
	 for( i=1 ; i< (nCellRow-1) ;++i)
	 {
	 	nInitAdh+=i*6;
	 }
	 
	 nInitAdh=floor(nInitAdh*nAdhPointCell*1.2);



/*
	for (i=0; i<nInitAdh ; ++i)
	{
		{
	 		List.push_back(Cadherin(lengthBal));
	 	//	InactiveCadh.push_back(i);
	 	}
	}
*/







/*
	 int i;

	 if (nCell==2)
	{
		for ( i=0; i<CadhSize; ++i)
		{
	 		List.push_back(Cadherin(2*i,2*i+1,lengthBal, initialCadherinNumber));
	 		ActiveCadh.push_back(i);
	 	}
	}

	else  // nCell==3
	{
		// Interface 1 and 2
		for ( i=0; i<CadhSize; ++i)
		{
	 		List.push_back(Cadherin(2*i,2*i+1,lengthBal, initialCadherinNumber));
	 		ActiveCadh.push_back(i);
	 	}

	 	// Cells 1 and 3
	 	int startD, startA;
	 	startD=(2*CadhSize-1)*3-4;
	 	startA=(2*CadhSize-1)*3-2;

	 	for ( i=0; i<CadhSize-1; ++i)
		{
	 		List.push_back(Cadherin(startD-4*i,startA-4*i,lengthBal, initialCadherinNumber));
	 		ActiveCadh.push_back(i);
	 	}

	 	i=CadhSize-1;
	 	List.push_back(Cadherin(startD-4*i+1,startA-4*i+1,lengthBal, initialCadherinNumber));
	 	ActiveCadh.push_back(i);

	 	// Cells 2 and 3
	 	startD=(2*CadhSize-1)*3-3;
	 	startA=(2*CadhSize-1)*3-1;

	 	for ( i=0; i<CadhSize-1; ++i)
		{
	 		List.push_back(Cadherin(startD-4*i,startA-4*i,lengthBal, initialCadherinNumber));
	 		ActiveCadh.push_back(i);
	 	}
	 	i=CadhSize-1;
	 	List.push_back(Cadherin(startD-4*i+1,startA-4*i,lengthBal, initialCadherinNumber));
	 	ActiveCadh.push_back(i);
	} 	

	*/

	// AdhCt.resize(List.size(),double(0));	
	 
  	//***********************************************************************// Update when 3 cells
 }
void CadherinSet::AddFreeCadherin()
{
	int ListID;

	List.push_back(Cadherin(lengthBal));
 	ListID=List.size()-1;


 	InactiveCadh.push_back(ListID);
 	AdhCt.push_back(0);

 	 List.at(ListID).firstAdhPoint=0;
	 List.at(ListID).secondAdhPoint=1;
	 List.at(ListID).force=double(0);
 	 List.at(ListID).ubProb=double(0);

	
	// List.at(ListID).gapCtAux=true;

 	List.at(ListID).active=false;
 	List.at(ListID).full=false;
	List.at(ListID).cadherinNumber=0;

}

 void CadherinSet::BindCadherin(int pointA, int pointB, int nInitAdh)
 {
 	int ListID;
 	//ListID=InactiveCadh[InactiveCadh.size()-1];

 	List.push_back(Cadherin(lengthBal));
 	ListID=List.size()-1;

 //	InactiveCadh.pop_back();
 	ActiveCadh.push_back(ListID);
 	AdhCt.push_back(0);

 	 List.at(ListID).firstAdhPoint=pointA;
	 List.at(ListID).secondAdhPoint=pointB;

 	List.at(ListID).active=true;
	List.at(ListID).cadherinNumber = nInitAdh ;

 }

 void CadherinSet::BindCadherinNew(int pointA, int pointB, int ListID, double distance)
 {
 	//int ListID;
 	//ListID=InactiveCadh[InactiveCadh.size()-1];

	List.at(ListID).length=distance;
 	ActiveCadh.push_back(ListID);
 	
 	List.at(ListID).firstAdhPoint=pointA;
	List.at(ListID).secondAdhPoint=pointB;

 	List.at(ListID).active=true;
	List.at(ListID).cadherinNumber = 1 ;


 }

 void CadherinSet::updateCts(int timeSlots)
 {
	for(int i=0; i<List.size(); ++i)
	{
		AdhCt.at(i)+=List.at(i).cadherinNumber/timeSlots;
	}
 }

void CadherinSet::UpdateCadhActivity()   // No usar ahora
{
	ActiveCadh.clear();
	InactiveCadh.clear();
	ActiveCadh.reserve(List.size());
	InactiveCadh.reserve(List.size());

	for(int i=0; i<List.size(); ++i)
	 {
	 	if ( List.at(i).active )
	 		ActiveCadh.push_back(i);	// Active is when it used for matrix
	 	if ( !List.at(i).full)
	 		InactiveCadh.push_back(i);  // Inactive is when is not full, only used for binding
	 }
}

void CadherinSet::UpdateUbProb(int ListID)
{
	double force;
	double K;
	force=List.at(ListID).force/List.at(ListID).cadherinNumber;

	K=KCatch*(exp(phiC-force/fAp)+exp(force/fAp-phiS));
	List.at(ListID).ubProb=1-exp(-K*TimeStep);
}

void CadherinSet::Binds(int ListID)  // When adhesion is full it does not get to this point
{

	List.at(ListID).active=true;
	if  (List.at(ListID).cadherinNumber == 0)
		ActiveCadh.push_back(ListID);
	else if  ( List.at(ListID).cadherinNumber ==(clusterMaxNumber-1) )
		List.at(ListID).full=true;
	else if ( List.at(ListID).cadherinNumber > (clusterMaxNumber-1) ) // If at the begginnig all of them are full, it is not stated therefore we need this to ensure that it works
	{
		List.at(ListID).full=true;
		List.at(ListID).cadherinNumber--; 
	}
	List.at(ListID).cadherinNumber++;
}

void CadherinSet::Unbinds(int ListID)
{
	if  (List.at(ListID).cadherinNumber == 1)
	 {	
	 	List.at(ListID).force=double(0);
 	 	List.at(ListID).ubProb=double(0);
	 	//List.at(ListID).length=lengthBal;
	 	List.at(ListID).active=false;
	 	List.at(ListID).gapCtAux=true;
	 	List.at(ListID).nGap++;
	 	List.at(ListID).gapGenTime.push_back(double(0));
	 	List.at(ListID).firstAdhPoint=0;
	 	List.at(ListID).secondAdhPoint=1;
	 }
	 List.at(ListID).full=false;
	 List.at(ListID).cadherinNumber--;
}

// *************** Cadh struct  ******************* //

Cadherin::Cadherin(int first, int second, int lengthCtr, int cadherinNumberInit)
 {
 	 force=double(0);
 	 ubProb=double(0);
	 length=lengthCtr;
	 firstAdhPoint=first;
	 secondAdhPoint=second;
	 active=true;
	 cadherinNumber=cadherinNumberInit;
	 nGap=0;
	 gapCtAux=false;
	 gapGenTime.push_back(double(0));
	 gapTotTime=double(0);
	 full=false;
 }

 Cadherin::Cadherin(int lengthCtr, int cadherinNumberInit)
 {
 	 force=double(0);
 	 ubProb=double(0);
	 length=lengthCtr;
	 firstAdhPoint=0;
	 secondAdhPoint=0;
	 active=false;
	 cadherinNumber=cadherinNumberInit;
	 nGap=0;
	 gapCtAux=false;
	 gapGenTime.push_back(double(0));
	 gapTotTime=double(0);
	 full=false;
 }

Cadherin::Cadherin(int lengthCtr)
{
 	 force=double(0);
 	 ubProb=double(0);
	 length=lengthCtr;
	 firstAdhPoint=0;
	 secondAdhPoint=0;
	 active=false;
	 cadherinNumber=0;
	 nGap=0;
	 gapCtAux=false;
	 gapGenTime.push_back(double(0));
	 gapTotTime=double(0);
	 full=false;
}