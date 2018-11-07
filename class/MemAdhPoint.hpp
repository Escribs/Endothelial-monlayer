// ##################################################
// #   MemAdhPoint.hpp - finally revised on Jun 2018     #
// #   coded by Jorge Escribano                     #
// #   Copyright (C) 2016-2018, Jorge Escribano,    #
// #   All rights reserved.                         #
// ##################################################

#ifndef MemAdhPoint_HPP
#define MemAdhPoint_HPP

#include "../header.hpp"
#include "../config/SingleConfig.hpp"

class MemAdhPoint
{
	public:
		MemAdhPoint(std::vector<double> CoordCtr, int IDCtr, int nAdhPerPoint);
	//	~ MemAdhPoint();
		
		int get_IDNode(){return IDNode;};
		int get_bound(){return bound;};
		bool GETFull(){return full;};
		bool GETGapIncluded(){return gapIncluded;};
		void SETBdry();

		void SETbound(bool a){bound=a;};
		void AddnCadhBound(int a){nCadhBound+=a;}
		
		void SETCadhBound(int IDCadh);
		void SETCadhFree(int IDCadh);
		
		void SETGapIncluded(bool a){gapIncluded=a;};

//		int GETcadhPointer(){return cadhPointer;};
//		void SETCadhBound(int a){cadhPointer=a;};

		std::vector<double> Coord;   	// Position of cells
		std::vector<int> CadhPointer;
		std::vector<int> ClosestPoints;

		// Plot variables
		int nGap;						// Number of times that a gap is generated
		bool gapCtAux;					// Only true at the step that it completely unbinds
		double gapTotTime;				// Total time withou adhesion
		std::vector<double> gapGenTime; // For each nGap



	private:
		bool bound;
		int IDNode;
		int nCadhBound;
		bool full;
		int nCadhMax;   // Maximun number of cadherins that can adhere to it
		bool gapIncluded;  // If is already contained in GAP
		// Poner 
};


#endif /* MemAdhPoint_HPP */