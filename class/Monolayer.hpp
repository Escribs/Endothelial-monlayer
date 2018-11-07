// ##################################################
// #   Monolayer.hpp - finally revised on Jun 2018     #
// #   coded by Jorge Escribano                     #
// #   Copyright (C) 2016-2018, Jorge Escribano,    #
// #   All rights reserved.                         #
// ##################################################

#ifndef Monolayer_HPP
#define Monolayer_HPP

#include "Cell.hpp"
#include "CadherinSet.hpp"

struct MembAdhData  // Structure to record adhesion charasteristic 
{
	MembAdhData();
	double distance;  // Current distance
	
	int refID;		// Number of the cell to which it is bound	
	
	int pointMinID;
	bool bound;						// If cadherin 

};

struct DistanceData  // Structure to record adhesion charasteristic 
{
	DistanceData();
	double distance;  // Current distance

	int firstPoint, secondPoint;
};

struct GapType
{
	GapType();

	int IDGapType;  // Number in the vector

	bool vertex;
	int stepsActive;
	int startStep;
	int endTime;
	bool active;
	int closestVert; // 0 if start points is closest or 1 if end point is closest

	int ctGapsID;




	double Area;
	int InitPoint;  // 
	int LastPoint;
	double Length;
	int nAdhPoint;
	double HeightAver;
	double HeightMax;
	double DistVertex;
	std::vector<int> IDList;   // List with the adh points. It does not have to ber ordered
	std::vector<int> CurrentShapeIDList;   // List of all the nodes that conforms the gap

};



class Monolayer
{
public:
	 Monolayer();
	//~ Monolayer();

	// Global equations
	 void CheckReBinding();
	 void Binding();
	 double BindingEqNew(int a, int b);     // If binding returns length if no -1
	 void Balance();
	 void printVTK(int ctParaview);
	 void UpdateStep(){step+=1;};
	 void printDataGlobal(char *index_cluster);
	 void printDataStep(char *index_cluster);
	 void matlabPlot(char *index_cluster);
	 void MakeUnboundForceZero(); // Optimize with unbound



	 void UpdateCadhProp();

	// Inside global, system equations
	 void ForceGeneration();
	 //SpMat BuildGlobalK();

	 void BuildSparseMatrixMembBdry();

	 //void TripletXYMemb(std::vector<T>& KTriplet,std::vector<T>& CTriplet,  int a, int b);
	// void TripletXYBdry(std::vector<T>& KTriplet,std::vector<T>& CTriplet,  int a, int b);
	// void TripletXYCadh(std::vector<T>& KTriplet,std::vector<T>& CTriplet,  int a, int b);

	 void TripletXY(std::vector<T>& KTriplet,std::vector<T>& CTriplet,  double EA, 
							double visc, double restLength, int a, int b, bool bdry);
	 void BuildTripletLocalKBdry( std::vector<T>& tripletList, int b, double coeff, double c, double s);
	 void BuildTripletLocalK( std::vector<T>& tripletList, int a, int b, double coeff, double c, double s);

	 void BuilSparseMatrixCadhConnectivity(SpMat &KMat, SpMat &CMat);
	 void BuildTripletLocalK( std::vector<T>& tripletList, int a, int b, double coeff, bool bdry);
	 void CheckUnbinding();
	 void CalculateForceIntern(Eigen::VectorXd &FInt);
	 void CalculateForceBdry(Eigen::VectorXd &FInt);
	 void CalculateForceMemb(Eigen::VectorXd &FInt);
	 void CalculateForceCadh(Eigen::VectorXd &FInt);
	 void CalculateForceBendingMemb(Eigen::VectorXd &FInt);
	 void CalculateForceBendingCadh(Eigen::VectorXd &FInt);
	 void CalculateNodesForce(Eigen::VectorXd &FInt, double EABalance, double barVisc, double lBalance , int a, int b);
	 void CalculateNodesForceBdry(Eigen::VectorXd &FInt, double EABalance, double barVisc, double lBalance , int a, int b);
	 void CalculateNodesForceCadh(Eigen::VectorXd &FInt, double EABalance, double barVisc, double lBalance , int a, int b); // Included to add repusion between cells

	 void RepForceCenter(Eigen::VectorXd &FInt);	// Force acting omn center by distance beetween centers

	 double CalculateStress(Eigen::VectorXd &FInt, double EABalance, double lBalance , int a, int b);
	 void RepulsionCellInCenter(Eigen::VectorXd &FInt);


	 void UpdateNodesPosition(Eigen::VectorXd &XIncEfeective);
	 void UpdateNodesPosition(std::vector<double> &XIncEfeective);
	 void CheckNodesIncrement();
	 void NodesInitVDir(int aIndex, int bIndex, std::vector<double> &VDir); // b-a
	 double KMatrixCoef(double lCurrent,double lBalance, double EABalance);

	 void CalculateForceInternAlternative(SpMat &KMat, Eigen::VectorXd &FInt, Eigen::VectorXd &XIncTotal );

	 void BalanceToy();
	 void print_data_intern(Eigen::MatrixXd &Amat, SpMat K_matrix_cc, Eigen::VectorXd &RNodes, Eigen::VectorXd &FInt , 
  								Eigen::VectorXd &Damp, Eigen::VectorXd &XIncTotal,	Eigen::VectorXd &XIncIt, int iter);
	 void BalanceLangevin();
	 void CalculateDamper(Eigen::VectorXd &CInt);
	 void CalculateNodesDamper(Eigen::VectorXd &CInt, double viscCoef , int a, int b, bool bdry);
	 void ApplyBdryConditions(Eigen::VectorXd &VectorAux);
	 void BindCadherin(int pointA, int pointB, int nInitAdh, int cadhID);
	 void PrintGAP(char *index_cluster);

	 void printDataGlobalVertex(char *index_cluster);
	 void printDataStepVertex(char *index_cluster, int t);

	 void SetInteractionsCells();
	 void SetInteractionsFace(int cellIDa, int cellIDb);

	 int GETVertexID(int cellID, int vertex);

	// General hepful equation
	void VDirCalculation(std::vector<double> a, std::vector<double> b, std::vector<double> &VDir);
	double Vmod(std::vector<double> a, std::vector<double> b);
	double NodesDist(int aIndex, int bIndex); // b-a
	double NodesDistX(int aIndex, int bIndex); // b-a
	double NodesDistY(int aIndex, int bIndex); // b-a
	double NodesDistInit(int aIndex, int bIndex);  // b-a
	double NodesDistInitX(int aIndex, int bIndex); // b-a
	double NodesDistInitY(int aIndex, int bIndex); // b-a
	void NodesVDir(int aIndex, int bIndex, std::vector<double> &VDir);
	void AngleVDir(double angle, std::vector<double> &VDir); 
	double NodesAngleInit(int aIndex, int bIndex);  // b-a
	double NodesAngle(int aIndex, int bIndex); // b-a
	void CalculateMiddlePoint(int aIndex, int bIndex, std::vector<double> &middlePoint); // b-a
	double AngleBalanceCenter(int a, int b);
	void CalculateAngleInit(Eigen::VectorXd &FInt);
	void UpdateMembpointFree(int aIndex,int bIndex, int IDCadh);
	void CheckUnboundAdhVectors();
	void InitialCheckofBoundSate();
	void AddSetFreeCadherins(int nCadhFree);
	void UpdateAdhPointCts();
	void CalculateMembraneDistance();
	double MinimumMembraneDistance(int refPointID, int cellID,int &pointMin);  // returns Distance between refPoint and closest membrane point of cellID
	int NormalizedID(int ID, int start); // Start is the first point of the cell (GetVertex(nCell, 0))
	void CalculateGapArea(GapType &GapAux);
	double polygonArea(std::vector<int>& pointsList); 

	bool IsPoligonSymple(std::vector<int> PointList );
	double AreaPoligonSymple(std::vector<int> PointList );

	void GetBoundListOK(std::vector<int> &ListAux , int ID );
	
	void CheckGapFormation(int step);
	void UpdateGaps(int step);
	double MinumunDistanceToVertex(std::vector<int> ListAux, int start, int end, int cellID, int& closestVertAux);
	void UpdateGAPVector();

	void PrintGAPStep(char *index_cluster);
	bool CheckGapOverlap(GapType &GapAux);
	void PrintGlobalStat(char *index_cluster);

	void ContractionByRepulsion(Eigen::VectorXd &FInt, int ID, double multiplier);
	double GapStress(GapType &GapAux );
	void SetBoundBdryConditions();




	void abort_program();


	void CalculateBendingForce(Eigen::VectorXd &FInt, double KBending, double barVisc, double angleBalance , int a, int b);
	void CalculateSFOverlapForce(Eigen::VectorXd &FInt , int centerP, int membP, int cellIDAux);
	void SFOverlapForceVTK(Eigen::VectorXd &FInt);

	void print_data_intern(SpMat K_matrix_cc, Eigen::VectorXd &RNodes, Eigen::VectorXd &FInt , 
  						Eigen::VectorXd &Damp, Eigen::VectorXd &XIncTotal,	Eigen::VectorXd &XIncIt,int divider, int iter);
	void print_data_intern( Eigen::VectorXd &FInt , Eigen::VectorXd &CInt, 
  									Eigen::VectorXd &XIncStep,	Eigen::VectorXd &Xvisc, int iter);

 	void print_data_intern2( Eigen::VectorXd &FNew );

 	void print_data_intern2(Eigen::VectorXd &FInt,Eigen::VectorXd &Faux , Eigen::VectorXd &Fvisc,
  						Eigen::VectorXd &XFext, Eigen::VectorXd &Xvisc, std::vector<double> &XIncIter, Eigen::VectorXd &CInt);

	void print_data_intern2(Eigen::VectorXd &FInt,Eigen::VectorXd &Faux , Eigen::VectorXd &Fvisc,
  						Eigen::VectorXd &XFext, Eigen::VectorXd &Xvisc, std::vector<double> &XIncIter, 
  						Eigen::VectorXd &XIncLastStep,int iter, int totalIter, Eigen::VectorXd &CInt);


	double CellOverlap(Eigen::VectorXd &FInt, int a, int b, double multiplier);

	int CellIDfromPoint(int a){return floor(a/(nAdhPointCell+1)); }
	bool IsMajorInCell(int a, int b); // Return true if a is major
	int ReturnMajorInCell(int a, int b); // Return major point


	void RepulsionMembraneGlobalForce(Eigen::VectorXd &FInt);
	void RepulsionMembraneGlobalInteraction(Eigen::VectorXd &FInt);
	void CheckMembDistanceGlobal(int IDaux, int IDcell, std::vector<bool> &MembDistGlobChecked);
	double CellClosestPoint(int IDref, int IDcellCheck, int &IDclosestAux);
	bool CheckOldGap(GapType &auxGap); // Return true if is equal to old gap. 
	void CellRemodel();

	void PrintStepStat(char *index_cluster);



	// Getters
	 int get_step(){return step;};
	 bool  GETprintFilesVTK(){return printFilesVTK;};
	 int GETTimeSlots(){return TimeSlots;};
	 int GETplotFreq(){return plotFreq;};

	// Variables
	 Eigen::VectorXd ForceLastStep;
	 Eigen::VectorXd ForceStep;
	 Eigen::VectorXd ForceInc;
	 Eigen::VectorXd ForceInternLastStep;
	 Eigen::VectorXd XIncTotal;
	 Eigen::VectorXd XIncTotalCheck;  // Measured in the nodes (node-node_init_pos)
	 Eigen::VectorXd XIncLastStep;
	 Eigen::VectorXd FDamperAcum;
	 Eigen::VectorXd CIntGlob;
	 Eigen::VectorXd ForceStepBr;

	 std::vector<MembAdhData> MembDistance;   // This vector corresponds to the minimum distance of the cell center mebrane points to other cell
	 									 	// Points ordered by ID

	 std::vector<DistanceData> MembDistanceGlobal;   // This vector corresponds to the minimum distance of the cell center mebrane points to other cell

	 									 	// Points ordered by ID
	 //std::vector<T> tripletListMembBdry; 
	 SpMat KMemBdry;
	 SpMat CMemBdry;
	std::vector<int> BdryPrint;      // To plot in vtk. 1 direction impeded; 0 -> free
	std::vector<std::vector<double>> NodesInitialPos;		// vecor with nodes initial position
	std::vector<std::vector<int>> ConnectivityMemb;		//Conenectivity to membrane points (first, second)
	std::vector<std::vector<int>> ConnectivityBdry;		//Connectivity to center point (center, memb)
	std::vector<Cell> Pool;								// Vector of cells
	std::vector<int> BdryIDs;// In global coordintas considering x and y (Apply Disp=0 )


	std::map<int, std::shared_ptr<MemAdhPoint>>	Nodes;
	std::map<int, std::shared_ptr<MemAdhPoint>>	NodesMemb;
	std::map<int, std::shared_ptr<MemAdhPoint>>	NodesBdry;
	std::vector<int> MoveableNodes;   // List of nodes that are not encastred

	std::vector<double> SFlength; 
	std::vector<double> SFlengthInit; 

	std::vector<double> MembLength; 
	std::vector<double> MembLengthInit; 


	std::default_random_engine generator; // For random time

//	std::vector<CadherinSet> CadhPack;   // One membEAr vecor, so I initialize it when I want

	std::vector<GapType> GAP;// Active Gap
	std::vector<GapType> GAPStepClosed;// Active Gap
	
	std::vector<int> InactiveGAP;  // Store inactive Gap
	std::vector<int> ActiveGAP;  // Store inactive Gap
	CadherinSet CadhPack;
	
	char *index_cluster_Aux;



	private:

	bool cellCenter,testFormulation;
	bool printFilesVTK;
	bool generalRupture;

	int ctVertexGlob;
	int ctInterfaceGlob;
	int ctVertexFirst;
	int ctInterfaceFirst;

	double percentGapDissap;
	double maximumRemodel;
	double centerMinimumDist;
	double cellAreaInit;   // Sum of stress fibers length

	int plotFreq;
	int meshSize;
	int nAdhPointSide;	// Number of adhesive points per side
	int nAdhPointCell;  // Number of adhesive points per cell
	int nCenterPointCell; 	// Number of center point is cell
	int stepForce;     // Number of time steps until force is recalculated
	int stepSmoothForce;
	int forceCt; 		// Counter for force aplication

	double GapMinTime;

	int nArists;
	int nCell;			// Total number of cells
	int nCellRow;
	double sepAdhPoint, sepInterface;
	double CadhLength;
	int step;
	double TimeStep;
	int TimeSlots;
	int smallDefResolution;
	int extraBdry;
	// For Gap track amd plot
	int ctVertexPrint;
	int ctInterfacePrint;
	int ctVertexStep;
	int ctInterfaceStep;

	int stepBrownian;


	double aristLength;
	double mediumViscCoef;
	double nucleusViscCoef;
	double membEA;
	double bdryEA;
	double membVisc;
	double bdryVisc;
	double maxCtForce;
	double maxRadForce;
	double maxFlowForce, dirFlowForce;
	double maximunDispAllowed;
	double KBendingMemb;
	double vertexLimitDistance;
	double gapThreshold;
	double maxForceRep;
	double cellRealArea;

	int nPropagationMemb, nPropagationBdry; 
	double ratioZeroRad, ratioZeroCt;  

	       // Remodeling

       double cellExtraGrowth, Kremodel, remodelRatioVel, remodelRatioVelVS;
       bool remodelUB;

       // Protrusion

        double maxProtForce, ratioZeroProt, stepForceProtr;
        int nPropagationProt;


};



#endif /* Monolayer_HPP */