// ##################################################
// #   MemAdhPoint.cpp - finally revised on Jun 2018       #
// #   coded by Jorge Escribano                     #
// #   Copyright (C) 2016-2018, Jorge Escribano,    #
// #   All rights reserved.                         #
// ##################################################

// Membrane points for adhesion class

#include "MemAdhPoint.hpp"


MemAdhPoint::MemAdhPoint(std::vector<double> CoordCtr, int IDCtr, int nAdhPerPoint)
 {
 	  IDNode=IDCtr;
  	std::vector<double>::iterator it;
  	it=CoordCtr.begin();
  	Coord.assign (it,CoordCtr.end());
    CadhPointer.resize(nAdhPerPoint,-1);
    nCadhMax=nAdhPerPoint;
    nCadhBound=0;
  	bound=false;
  	full=false;

  			// Plot variables
	 nGap=0;						// Number of times that a gap is generated.
	 gapCtAux=false;					// Only true at the step that it completely unbinds
	 gapGenTime.push_back(double(0));
	 gapTotTime=double(0);
   gapIncluded=false;

  	
 }

  void MemAdhPoint::SETBdry()
  {

    bound=true;
    full=true;

  }

 void MemAdhPoint::SETCadhBound(int IDCadh)
 {
 	nCadhBound++;
 	CadhPointer.at(nCadhBound-1)=IDCadh;
 	bound=true;
 	if (nCadhMax==nCadhBound)
 		full=true;
 	else
 		full=false;
 }

 void MemAdhPoint::SETCadhFree(int IDCadh)
 {
 	nCadhBound--;
 	full=false;
 	int auxID;
 	if (nCadhBound > 0)
 		bound=true;
 	else
 	{
 		bound=false;

 		gapCtAux=true;
	 	nGap++;
	 	gapGenTime.push_back(double(0));
 	}


 	for (auto it= CadhPointer.begin(); it!=CadhPointer.end();)
 	{
 		auxID=*it;
 		if (auxID== IDCadh)
 		{
 			CadhPointer.erase(it);
 			it=CadhPointer.end();
 			CadhPointer.push_back(-1);
 		}

 		else
 		it++;

 	}
 }