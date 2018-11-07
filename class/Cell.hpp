// ##################################################
// #   Cell.hpp - finally revised on Jun 2018     #
// #   coded by Jorge Escribano                     #
// #   Copyright (C) 2016-2018, Jorge Escribano,    #
// #   All rights reserved.                         #
// ##################################################

#ifndef Cell_HPP
#define Cell_HPP

#include "MemAdhPoint.hpp"

class Cell
{
	 public:
	 	Cell( std::vector<double> Center ,std::vector<std::vector<int>>& ConnectivityMemb, std::vector<std::vector<int>>& ConnectivityBdry, int IDCtr);
	 	Cell(std::vector<std::vector<int>>& ConnectivityMemb, std::vector<std::vector<int>>& ConnectivityBdry, int IDCtr);
	 //	Cell(int IDCtr, bool cellCenter, int nArists);
    	void PointsDirection(std::vector<double> &Coord, std::vector<double> &VDir, bool Membrane);
    	int GETnArists(){return nArists;};
    	int GETnAdhPoint(){return nAdhPoint;};
	    int GETIDCenter(){return IDCenter;};

    	void PointsDirection3cells( std::vector<double> &VDir);
		void POINTSCHECK();
	 	
		std::vector<std::shared_ptr<MemAdhPoint>> PointMemb;		// Membrane adhesive points
		std::vector<std::shared_ptr<MemAdhPoint>> PointBdry;		// Boundary condition points

		std::set<int> UnboundAdhPoint;		// Vector containing Point memb ub

		std::vector<double> CenterPoint;



	 private:
		bool cellCenter; // Determines boundary condition; 1 -> Center condition | 2 -> Wall condition

		int ID; 	// Difference each cell 1 left , 2 right, 3 top
		int nArists;
		int nAdhPoint;		//Number of adhesive points
		int nAdhPointSide;		//Number of adhesive points
		int IDCenter;

		double aristLength;
		double sepAdhPoint; // Separation between adhesive points
		double sepInterface;
		double cadhLength;
		double radius;
}; 

#endif /* Cell_HPP */