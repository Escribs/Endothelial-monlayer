// ##################################################
// #   Monolayer.cpp - finally revised on Jun 2018  #
// #   coded by Jorge Escribano                     #
// #   Copyright (C) 2016-2018, Jorge Escribano,    #
// #   All rights reserved.                         #
// ##################################################

// Main class

#include "Monolayer.hpp"
MembAdhData::MembAdhData()
{
	 distance=0;  // Current distance
	
	 refID=0;		// Number of the cell to which it is bound	

	pointMinID=0;
	 bound=false;						// If cadherin 
}

DistanceData::DistanceData()
{
 	distance=double(0);  // Current distance
	firstPoint=0;
	secondPoint=0;
}


GapType::GapType()
{
 
	 stepsActive=0;
	 startStep=0;
	 endTime=0;
	 closestVert=0;
	 Area=0.;
	 InitPoint=0;  // 
	 LastPoint=0;
	 Length=0.;
	 nAdhPoint=0;
	 HeightAver=0.;
	 HeightMax=0.;
	 DistVertex=0.;

}



 Monolayer::Monolayer()
 {
 	int i,j,z, cellCt;
 	double xlengthAux, ylengthAux;
 	generator.seed( std::chrono::high_resolution_clock::now().time_since_epoch().count() );
 	std::vector<double> Coord;
 	std::vector<double> CoordRef;
 	std::vector<int> BdryAux;

	// Read variables
	 SingleConfig DataConfig=SingleConfig::getConfig();

	 nucleusViscCoef=DataConfig.cadhVisc;   // Uso la de la cadherina por la cadhesrina esta comentada y asi no modifico input
	 										// Mirar en CalculateDamper


	 nCellRow=DataConfig.nCell;
	 TimeStep=DataConfig.TimeStep;
	 TimeSlots=DataConfig.TimeSlots;
	 membEA=DataConfig.membEA;
	 bdryEA=DataConfig.bdryEA;
	 membVisc=DataConfig.membVisc;
	 bdryVisc=DataConfig.bdryVisc;
	 sepAdhPoint=DataConfig.sepAdhPoint;
	 sepInterface=DataConfig.sepInterface;
	 maxCtForce=DataConfig.maxCtForce;
 	 nArists=DataConfig.nArists;
 	 cellCenter=DataConfig.cellCenter;
 	 testFormulation=DataConfig.testFormulation;
 	 smallDefResolution=DataConfig.smallDefResolution;
 	 mediumViscCoef=DataConfig.mediumViscCoef;
 	 aristLength=DataConfig.aristLength;
 	 printFilesVTK=DataConfig.printFilesVTK;
 	 plotFreq=DataConfig.plotFreq;

 	 stepBrownian=0;



 	 maxRadForce=DataConfig.maxRadForce;
 	 maxFlowForce=DataConfig.maxFlowForce;
 	 dirFlowForce=DataConfig.dirFlowForce;
 	 extraBdry=DataConfig.extraBdry;
 	 KBendingMemb=DataConfig.KBendingMemb;


 	 nPropagationMemb=DataConfig.nPropagationMemb; 
 	 nPropagationBdry=DataConfig.nPropagationBdry; 
	 ratioZeroRad=DataConfig.ratioZeroRad;
	 ratioZeroCt=DataConfig.ratioZeroCt;  

 	 maximunDispAllowed=DataConfig.maximunDispAllowed;

 	 gapThreshold=DataConfig.gapThreshold*sepInterface*sepAdhPoint;
 	 int initialCadherinNumber=DataConfig.initialCadherinNumber;

 	 
 	 remodelRatioVel=DataConfig.remodelRatioVel;
 	 Kremodel=DataConfig.Kremodel;
 	 cellExtraGrowth=DataConfig.cellExtraGrowth;
 	 remodelRatioVelVS=DataConfig.remodelRatioVelVS;
 	 remodelUB=DataConfig.remodelUB;


 	 maxProtForce=DataConfig.maxProtForce;
 	 ratioZeroProt=DataConfig.ratioZeroProt;
 	 stepForceProtr=DataConfig.stepForceProtr;
 	 nPropagationProt=DataConfig.nPropagationProt;
 
 	 maxForceRep=maxRadForce+maxCtForce;
 	 if ( maxForceRep < 0.04)  
 	 	maxForceRep=0.04;


 	 percentGapDissap=0.7;

 	 maximumRemodel=(maxForceRep)/bdryVisc*TimeStep*Kremodel*remodelRatioVel;





	 ctVertexPrint=0;
	 ctInterfacePrint=0;
	 ctVertexStep=0;
	 ctInterfaceStep=0;
	 ctInterfaceFirst=0;
	 ctVertexFirst=0;
 	
 	 stepSmoothForce=DataConfig.stepSmoothForce;
 	 stepForce=DataConfig.stepForce;
	 forceCt=stepForce+1;  // In the firs step enters the loop to claculate force
 	 nArists=6;
	 step=0;

	// Create Cells
	 ConnectivityBdry.resize(2);
	 ConnectivityMemb.resize(2);
	 std::vector<std::shared_ptr<MemAdhPoint>>::iterator it;

	// Center points

	 Coord.push_back(double(0));
	 Coord.push_back(double(0));
	 CoordRef=Coord;
	 nCenterPointCell=1;

	 nAdhPointSide=floor((aristLength+0.1)/sepAdhPoint); 
	 aristLength=double(nAdhPointSide*sepAdhPoint);
	 nAdhPointCell=double(nAdhPointSide*nArists);
	 MembDistance.resize(nAdhPointCell);



	 // Distance between centers
	 xlengthAux= 2*aristLength*cos(PI/6)+sepInterface;
	 ylengthAux= aristLength/2+sepInterface/2*tan(PI/3)+aristLength;

	 centerMinimumDist=2*aristLength*cos(PI/6);  // For repuslion between centers

	 vertexLimitDistance=aristLength/5;
	 generalRupture=false;


	 int row, column;

	 column=nCellRow;
	 row=nCellRow;
	 cellCt=0;
	 
	 double xAux;
	 // Center and upper zone
	 for (j=0; j < column; ++j)
	 {
	 	row=nCellRow-j;
	 	xAux=CoordRef[0]+xlengthAux/2*j;
	 	Coord[1]=CoordRef[1]+j*ylengthAux;
	 	for (i=0; i < row; ++i)
	 	{
			Coord[0]=CoordRef[0]+xAux+xlengthAux*i;
		 	Pool.push_back(Cell(Coord,ConnectivityMemb, ConnectivityBdry, cellCt));

		 	cellCt++;
	 	}	 
	}

	// Lower zone
	 row=nCellRow;
	 for (j=1; j < column; ++j)
	 {
	 	row=nCellRow-j;
	 	xAux=CoordRef[0]+xlengthAux/2*j;
	 	Coord[1]=CoordRef[1]-j*ylengthAux;
	 	for (i=0; i < row; ++i)
	 	{
			Coord[0]=xAux+CoordRef[0]+xlengthAux*i;
		 	Pool.push_back(Cell(Coord,ConnectivityMemb, ConnectivityBdry, cellCt));
		 	cellCt++;
	 	}	 
	 }
	
	  nCell=cellCt;
	// Copy points to nodes  
	 for (i=0;i< Pool.size();++i)	
	 {	 	//for (int j=0;Pool.at(i).PointMemb.size(); ++j)	
		 		//Nodes[Pool.at(i).PointMemb.at(j)->get_IDNode()]=Pool.at(i).PointMemb.at(j);
		for (it=Pool.at(i).PointMemb.begin();it != Pool.at(i).PointMemb.end();++it)
		 {
		 	Nodes[(*it)->get_IDNode()]=(*it);
		 	NodesMemb[(*it)->get_IDNode()]=(*it);
		 }
		 for (it=Pool.at(i).PointBdry.begin();it != Pool.at(i).PointBdry.end();++it)
		 {
		 	Nodes[(*it)->get_IDNode()]=(*it);
		 	NodesBdry[(*it)->get_IDNode()]=(*it);
		 }
	 }

	// Apply boundary conditions
	 
	 cellCt=0;
	 row=nCellRow;
	 int encastredID, encastredIDmax;
	 // Center and upper zone
	 for (j=0; j < column-1; ++j)
	 {
	 	row=nCellRow-j;
	 	
	 	// bdry in x direction
		encastredID=GETVertexID(cellCt,4);
		encastredIDmax=GETVertexID(cellCt,5);
		
		for (i=encastredID; i<= encastredIDmax;++i )
		{
			BdryIDs.push_back(i*2);
			BdryIDs.push_back(i*2+1);
			BdryAux.push_back(i);
		
		}
		if (j==0)
		{
			// bdry in x direction
			encastredID=GETVertexID(cellCt,3);
			encastredIDmax=GETVertexID(cellCt,4);
		
			for (i=encastredID; i< encastredIDmax;++i )
			{
				BdryIDs.push_back(i*2);
				BdryIDs.push_back(i*2+1);
				BdryAux.push_back(i);

			}
		}

		// bdry in x direction
		encastredID=GETVertexID(cellCt,5);
		encastredIDmax=nAdhPointSide+encastredID;
		
		for (i=encastredID; i< encastredIDmax;++i )
		{
			BdryIDs.push_back(i*2);
			BdryIDs.push_back(i*2+1);
			BdryAux.push_back(i);

		}


		cellCt=cellCt+row;
		
		encastredID=GETVertexID(cellCt-1,0);
		encastredIDmax=GETVertexID(cellCt-1,1);

		for (i=encastredID; i< encastredIDmax;++i )
		{
			BdryIDs.push_back(i*2);
			BdryIDs.push_back(i*2+1);
			BdryAux.push_back(i);

		}
		

		// bdry in x direction
		encastredID=GETVertexID(cellCt-1,1);
		encastredIDmax=GETVertexID(cellCt-1,2);
		
		for (i=encastredID; i<= encastredIDmax;++i )
		{
			BdryIDs.push_back(i*2);
			BdryIDs.push_back(i*2+1);
			BdryAux.push_back(i);
		}

		if (j==0)
		{
			// bdry in x direction
			encastredID=GETVertexID(cellCt-1,2)+1;
			encastredIDmax=GETVertexID(cellCt-1,3);
		
			for (i=encastredID; i<= encastredIDmax;++i )
			{
				BdryIDs.push_back(i*2);
				BdryIDs.push_back(i*2+1);
				BdryAux.push_back(i);
			}
		}
	 }

	 // Upper cell
	 encastredID=GETVertexID(cellCt,4);
	 encastredIDmax=GETVertexID(cellCt,5);
		
	 for (i=encastredID; i<= encastredIDmax;++i )
	{
		BdryIDs.push_back(i*2);
		BdryIDs.push_back(i*2+1);
		BdryAux.push_back(i);

	}

	 
	 encastredID=GETVertexID(cellCt,1);
	 encastredIDmax=GETVertexID(cellCt,2);

	 for (i=encastredID; i<= encastredIDmax;++i )
	{
		BdryIDs.push_back(i*2);
		BdryIDs.push_back(i*2+1);
		BdryAux.push_back(i);
	}

	 encastredID=GETVertexID(cellCt,0);
	 encastredIDmax=GETVertexID(cellCt,1);

	 for (i=encastredID; i<= encastredIDmax;++i )
	{
		BdryIDs.push_back(i*2);
		BdryIDs.push_back(i*2+1);
		BdryAux.push_back(i);
	}

	 encastredID=GETVertexID(cellCt,5);
	 encastredIDmax=encastredID+nAdhPointSide;

	 for (i=encastredID; i< encastredIDmax;++i )
	{
		BdryIDs.push_back(i*2);
		BdryIDs.push_back(i*2+1);
		BdryAux.push_back(i);
	}



	// bdry in y direction
	 encastredID=GETVertexID(cellCt,0);
	 BdryIDs.push_back(encastredID*2+1);
	 cellCt++;

	 // Lower zone
	 row=nCellRow;
	 for (j=1; j < column-1; ++j)   // Check if it -1
	 {
	 	row=nCellRow-j;
		// bdry in x direction
		encastredID=GETVertexID(cellCt,4);
		encastredIDmax=GETVertexID(cellCt,5);
		
		for (i=encastredID; i<= encastredIDmax;++i )
		{
			BdryIDs.push_back(i*2);
			BdryIDs.push_back(i*2+1);
			BdryAux.push_back(i);

		}

		encastredID=GETVertexID(cellCt,3);
		encastredIDmax=GETVertexID(cellCt,4);
		
		for (i=encastredID; i< encastredIDmax;++i )
		{
			BdryIDs.push_back(i*2);
			BdryIDs.push_back(i*2+1);
			BdryAux.push_back(i);

		}


		cellCt=cellCt+row;

		encastredID=GETVertexID(cellCt-1,1);
		encastredIDmax=GETVertexID(cellCt-1,2);

		for (i=encastredID; i< encastredIDmax;++i )
		{
			BdryIDs.push_back(i*2);
			BdryIDs.push_back(i*2+1);
			BdryAux.push_back(i);
		}	
		encastredID=GETVertexID(cellCt-1,2);
		encastredIDmax=GETVertexID(cellCt-1,3);

		for (i=encastredID; i<= encastredIDmax;++i )
		{
			BdryIDs.push_back(i*2);
			BdryIDs.push_back(i*2+1);
			BdryAux.push_back(i);
		}	
	 }

	 // Lower cell
	 encastredID=GETVertexID(cellCt,4);
	 encastredIDmax=GETVertexID(cellCt,5);
		
	 for (i=encastredID; i<= encastredIDmax;++i )
	{
		BdryIDs.push_back(i*2);
		BdryIDs.push_back(i*2+1);
		BdryAux.push_back(i);
	}	
	 
	 encastredID=GETVertexID(cellCt,1);
	 encastredIDmax=GETVertexID(cellCt,2);

	 for (i=encastredID; i<= encastredIDmax;++i )
	{
		BdryIDs.push_back(i*2);
		BdryIDs.push_back(i*2+1);
		BdryAux.push_back(i);
	}	

	encastredID=GETVertexID(cellCt,2);
	 encastredIDmax=GETVertexID(cellCt,3);

	 for (i=encastredID; i< encastredIDmax;++i )
	{
		BdryIDs.push_back(i*2);
		BdryIDs.push_back(i*2+1);
		BdryAux.push_back(i);
	}	
	 	
	encastredID=GETVertexID(cellCt,3);
	 encastredIDmax=GETVertexID(cellCt,4);

	 for (i=encastredID; i< encastredIDmax;++i )
	{
		BdryIDs.push_back(i*2);
		BdryIDs.push_back(i*2+1);
		BdryAux.push_back(i);
	}	 	



	// Set initial cadherins 

	 int cadhCt=0;
	 int cellCtAux=0;
	 for (int s=1; s<= nCellRow;++s)
		cellCtAux+=s;


	 row=nCellRow;
	 cellCt=0;
	 int startA, startB;  // A it is always the lower in (y direction) cell
	 // Center and upper zone
	 for (j=0; j < column-1; ++j)
	 {
	 	row=nCellRow-j;
	 	
	 	for (i=0; i < row-1; ++i)
	 	{
	 		// horizontal position
	 		startA = GETVertexID(cellCt,1);
	 		startB = GETVertexID(cellCt+1,5);
	 		z=0;
	 		BindCadherin(startA+z, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
	 		cadhCt++;
			for (int z=1; z <nAdhPointSide ; z++ )
			{
				BindCadherin(startA+z, startB-z, initialCadherinNumber, cadhCt);
		 		cadhCt++;
			}
			z=nAdhPointSide;
			BindCadherin(startA+z, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
	 		cadhCt++;


			// Crossed in first angle postion

			startA = GETVertexID(cellCt,0);
	 		startB = GETVertexID(cellCt+row,4);
	 		z=0;
	 		BindCadherin(startA+z, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
	 		cadhCt++;
			for (int z=1; z <nAdhPointSide ; z++ )
			{
				BindCadherin(startA+z, startB-z, initialCadherinNumber, cadhCt);
	 			cadhCt++;
			}
			z=nAdhPointSide;
			BindCadherin(startA+z, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
	 		cadhCt++;

			// Crossed in second angle postion
			startA = GETVertexID(cellCt+1,5); // Change to 0 in last position
	 		startB = GETVertexID(cellCt+row,3);
	 		z=0;
	 		BindCadherin(startA+z, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
	 		cadhCt++;
			for (int z=1; z <nAdhPointSide ; z++ )
			{
				BindCadherin(startA+z, startB-z, initialCadherinNumber, cadhCt);
	 			cadhCt++;
			}
			z=nAdhPointSide;
			startA = GETVertexID(cellCt+1,0);
			BindCadherin(startA, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
	 		cadhCt++;

			if (j==0)
			{
				// Crossed in third angle postion

				
	
				startA = GETVertexID(cellCtAux+i,5);
		 		startB = GETVertexID(cellCt,3);
		 		z=0;
		 		BindCadherin(startA+z, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
	 			cadhCt++;
				for (int z=1; z <nAdhPointSide ; z++ )
				{
					BindCadherin(startA+z, startB-z, initialCadherinNumber, cadhCt);
			 		cadhCt++;
				}
				z=nAdhPointSide;
				startA = GETVertexID(cellCtAux+i,0);
				BindCadherin(startA, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
	 			cadhCt++;
	
				// Crossed in the fourth angle postion
	
				startA = GETVertexID(cellCtAux+i,0);
		 		startB = GETVertexID(cellCt+1,4);
		 		z=0;
		 		BindCadherin(startA+z, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
		 		cadhCt++;
				for (int z=1; z <nAdhPointSide ; z++ )
				{
					BindCadherin(startA+z, startB-z, initialCadherinNumber, cadhCt);
			 		cadhCt++;
				}
				z=nAdhPointSide;

				BindCadherin(startA+z, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
		 		cadhCt++;

			}

			cellCt++;
	 	}	
	 	cellCt++; 
	 }

	 cellCt++;

	 // Lower zone
	 row=nCellRow;
	 for (j=1; j < column-1; ++j)
	 {
	 	row=nCellRow-j;
	 	
	 	for (i=0; i < row-1; ++i)
	 	{
	 		// horizontal position
	 		startA = GETVertexID(cellCt,1);
	 		startB = GETVertexID(cellCt+1,5);
	 		z=0;
	 		BindCadherin(startA+z, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
	 		cadhCt++;
			for (int z=1; z <nAdhPointSide ; z++ )
			{
				BindCadherin(startA+z, startB-z, initialCadherinNumber, cadhCt);
	 			cadhCt++;
			}
			z=nAdhPointSide;
			BindCadherin(startA+z, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
	 		cadhCt++;


			// Crossed in fourth angle postion

			startA = GETVertexID(cellCt+row,5);
	 		startB = GETVertexID(cellCt,3);
	 		z=0;
	 		BindCadherin(startA+z, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
	 		cadhCt++;
			for (int z=1; z <nAdhPointSide ; z++ )
			{
				BindCadherin(startA+z, startB-z, initialCadherinNumber, cadhCt);
	 			cadhCt++;
			}
			z=nAdhPointSide;
			startA = GETVertexID(cellCt+row,0);
			BindCadherin(startA, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
	 		cadhCt++;

			// Crossed in second angle postion
			startA = GETVertexID(cellCt+row,0); // Change to 0 in last position
	 		startB = GETVertexID(cellCt+1,4);
	 		z=0;
	 		BindCadherin(startA+z, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
	 		cadhCt++;
			for (int z=1; z <nAdhPointSide ; z++ )
			{
				BindCadherin(startA+z, startB-z, initialCadherinNumber, cadhCt);
	 			cadhCt++;
			}
			z=nAdhPointSide;
			BindCadherin(startA+z, startB-z, floor(initialCadherinNumber/2)+1, cadhCt );
	 		cadhCt++;


			cellCt++;
	 	}	
	 	cellCt++; 
	 }



	 
	 //Save nodes initial position

	NodesInitialPos.resize(Nodes.size());
	i=0;
	for (auto itsec=Nodes.begin();itsec != Nodes.end();++itsec)
	{
		NodesInitialPos.at(i).resize(2);
		NodesInitialPos.at(i).at(0) =(itsec->second)->Coord.at(0);
		NodesInitialPos.at(i).at(1) =(itsec->second)->Coord.at(1);
		i++;
	}

	// Mesh initialazation, NodesMemb becuase boundry is fixed
	meshSize=Nodes.size();

	// Set bdry vector
	 BdryPrint.resize(2*meshSize, 0);

	 for (j=0; j<BdryIDs.size(); ++j)
	 	BdryPrint.at(BdryIDs[j])=1;

	 
//	BuildSparseMatrixMembBdry();

	ForceLastStep= Eigen::VectorXd::Zero(2*meshSize);
	ForceStepBr= Eigen::VectorXd::Zero(2*meshSize);
	ForceInc= Eigen::VectorXd::Zero(2*meshSize);
	XIncTotal= Eigen::VectorXd::Zero(2*meshSize);
	ForceStep= Eigen::VectorXd::Zero(2*meshSize);
//	ForceInternLastStep= Eigen::VectorXd::Zero(2*meshSize);
	XIncLastStep= Eigen::VectorXd::Zero(2*meshSize);
	FDamperAcum= Eigen::VectorXd::Zero(2*meshSize);
	XIncTotalCheck= Eigen::VectorXd::Zero(2*meshSize);


	InitialCheckofBoundSate();

	// Add free cadherins in case mora than one adhesion per node is included
	if (DataConfig.nAdhPerPoint > 1)
		AddSetFreeCadherins(cadhCt);


	CIntGlob= Eigen::VectorXd::Zero(2*meshSize);
	CalculateDamper(CIntGlob);

	MoveableNodes.reserve(meshSize);
	for (j=0; j<meshSize; ++j)
		MoveableNodes.push_back(j);

	// quitar los boundary

	int auxEncast;

	for (i=0; i<BdryAux.size(); ++i)
	{
		for (auto it=MoveableNodes.begin(); it!=MoveableNodes.end();)
		{
			auxEncast=*it;
			if (auxEncast==BdryAux[i])
			{
				MoveableNodes.erase(it);
				it=MoveableNodes.end();
			}
			else
				it++;	
		}

	}


	MembDistanceGlobal.reserve(nAdhPointCell*Pool.size());


	// Modify conecctivity membrane  
	SFlength.resize(Nodes.size(), double(0));
	MembLength.resize(Nodes.size(), double(0));
	for (i=0; i< ConnectivityBdry.at(0).size();++i )
	{
		
		SFlength.at(ConnectivityBdry.at(1).at(i))= NodesDistInit( ConnectivityBdry.at(0).at(i),ConnectivityBdry.at(1).at(i) );		

	}

	for (i=0; i< ConnectivityMemb.at(0).size();++i )
	{
		
		MembLength.at(ConnectivityMemb.at(0).at(i))= NodesDistInit( ConnectivityMemb.at(0).at(i),ConnectivityMemb.at(1).at(i) );		

	}
	

	SFlengthInit=SFlength;
	
	cellAreaInit=double(0);
	for (i=0;i < nAdhPointCell+1; ++i)
	{
		cellAreaInit+=SFlength[i];
	}

	cellAreaInit=cellAreaInit*cellExtraGrowth;

	cellRealArea= 2.6*aristLength*aristLength;

	SetBoundBdryConditions();
 }


void Monolayer::BindCadherin(int pointA, int pointB, int nInitAdh, int cadhID)
{
	CadhPack.BindCadherin( pointA, pointB, nInitAdh);
	
	// Make the AdhPoint part
	Nodes.at(pointA)-> SETCadhBound(cadhID);
	Nodes.at(pointB)-> SETCadhBound(cadhID);
	//Nodes.at(pointA)-> cadhPointer.push_back(cadhID);
	//Nodes.at(pointB)-> cadhPointer.push_back(cadhID);


}
void Monolayer::AddSetFreeCadherins(int nCadhFree)
{
	for (int i=0; i < nCadhFree; ++i)
		CadhPack.AddFreeCadherin();

}


void Monolayer::SetInteractionsCells()
{
	double celldistance;
    double distanceLimitCell=aristLength*3;

	// Cell center
	for (int i=0; i < nCell-1; ++i)
	{
		for (int j=i+1; j < nCell; ++j)
		{
			celldistance=NodesDist(Pool[i].GETIDCenter(), Pool[j].GETIDCenter());
			if (celldistance < distanceLimitCell )
				SetInteractionsFace(i,j);
		}
	}
}

void Monolayer::SetInteractionsFace(int cellIDa, int cellIDb)  
{
	
	double angleVal, distX, distY;
	int sideIDa, sideIDb;
	distX=NodesDistX(Pool[cellIDa].GETIDCenter(), Pool[cellIDb].GETIDCenter()); // b- a
	distY=NodesDistY(Pool[cellIDa].GETIDCenter(), Pool[cellIDb].GETIDCenter());
	angleVal=atan2(distY,distX);
	if ( angleVal < 0 ) 
		angleVal=angleVal+2*PI;

	if ( (angleVal > PI/6.) && (angleVal <= PI/2.) )  // Side 0
	{
		sideIDa=0;
		sideIDb=4;
	}
	else if ( (angleVal > PI/2.) && (angleVal <= 5*PI/6.) )  // Side 0
	{
		sideIDa=5;
		sideIDb=3;
	}	
	else if ( (angleVal > 5*PI/6.) && (angleVal <= 7*PI/6.) )  // Side 0
	{
		sideIDa=4;
		sideIDb=2;
	}	
	else if ( (angleVal > 7*PI/6.) && (angleVal <= 9*PI/6.) )  // Side 0
	{
		sideIDa=3;
		sideIDb=1;
	}
	else if ( (angleVal > 9*PI/6.) && (angleVal <= 11*PI/6.) )  // Side 0
	{
		sideIDa=2;
		sideIDb=0;
	}
	else
	{
		sideIDa=1;
		sideIDb=5;
	}		

}



void Monolayer::Binding()    // Check whole cell with whole cell
{
	double celldistance, adhPointDistance,adhPointDistanceAux;
    double distanceLimitCell=aristLength*3;
    int IDa, IDb;
    bool loopAux,checkBind;

	// Cell center
	for (int i=0; i < nCell-1; ++i)
	{
		for (int j=i+1; j < nCell; ++j)
		{
			celldistance=NodesDist(Pool[i].GETIDCenter(), Pool[j].GETIDCenter());
			if (celldistance < distanceLimitCell )
			{
				
				for (auto ita=Pool[i].UnboundAdhPoint.begin(); ita!=Pool[i].UnboundAdhPoint.end();)
				{
					IDa=*ita;
					adhPointDistance=-1;
					for (auto itb=Pool[j].UnboundAdhPoint.begin(); itb!=Pool[j].UnboundAdhPoint.end();)
					{
						IDb=*itb;
						loopAux=true;
						checkBind=false;
						// Funcion para chequear no repticion considerar si ponerla mas tarde
						//AvoidRepCadh()*******************************************************************************************************************

						// To check only the closest point
						adhPointDistance=NodesDist(IDa, IDb);
						while ( ( (adhPointDistance < CadhPack.GETbindLimit()) && (CadhPack.InactiveCadh.size()>0) ) &&   loopAux )
	 					{
	 						checkBind=true;
	 						loopAux=false;
	 						itb++;
	 						if ( itb != Pool[j].UnboundAdhPoint.end() )
	 						{
	 							IDb=*itb;
	 							adhPointDistanceAux=NodesDist(IDa, IDb);
	 							if (adhPointDistanceAux < adhPointDistance )
	 							{
	 								adhPointDistance=adhPointDistanceAux;
	 								loopAux=true;
	 							}
	 							else
	 							{
	 								itb--;
	 								IDb=*itb;
	 							}

	 						}
	 						else
	 						{
	 							itb--;
	 							IDb=*itb;
	 						}

	 					}																																																																																																									

	 					if (checkBind)
							adhPointDistance=BindingEqNew(IDa,IDb);
						else
							adhPointDistance=-1;

						if (adhPointDistance > 0 ) // Binding   They are bound
						{
							// Only erase of vector if the cadherin is full (update previously in BindingEqNew)
							if  ( Nodes.at(IDa) -> GETFull() )
								ita=Pool[i].UnboundAdhPoint.erase(ita);
							if  ( Nodes.at(IDb) -> GETFull() )
								Pool[j].UnboundAdhPoint.erase(itb);
							

							itb=Pool[j].UnboundAdhPoint.end();// We go to the next monomer

						}
						else
							itb++;

							

					}
					if (adhPointDistance < 0 )
						ita++;

				}
			}	
		}
	}


	CheckUnboundAdhVectors();    // To exclude from Ubound adhesion Point vector en each cell the nodes bound  
	  // Lo podemos poner si hay algun evento de union nuevo, si no quitarlo
}

double Monolayer::BindingEqNew(int a, int b)     // If binding returns length if no -1
{
	double adhPointDistance;
	std::uniform_real_distribution<double> distribution(0, 1);
	double prob, Kon, rnd;
	int IDCadh;

	adhPointDistance=NodesDist(a, b);
	if ( (adhPointDistance < CadhPack.GETbindLimit()) && (CadhPack.InactiveCadh.size()>0) )
	 {
	 	Kon=CadhPack.GETKonRate()*CadhPack.GETdensity()*
	 		(1-adhPointDistance/CadhPack.GETbindLimit());
	 	prob=(1-exp(-Kon*TimeStep));
	 	rnd=distribution(generator);
	 	if (prob > rnd) 
	 	{
	 		IDCadh=CadhPack.InactiveCadh[CadhPack.InactiveCadh.size()-1];
	 		CadhPack.InactiveCadh.pop_back();
	 		CadhPack.BindCadherinNew(a,b,IDCadh,adhPointDistance );
	 		Nodes.at(a) -> SETCadhBound(IDCadh);
			Nodes.at(b) -> SETCadhBound(IDCadh);
	 	}
	 	else
	 		adhPointDistance=-1;

	 }
	 else 
	 	adhPointDistance=-1;

	return adhPointDistance; 
}


void Monolayer::CheckReBinding()    //  Proportional to Distance
{
	std::uniform_real_distribution<double> distribution(0, 1);
	double prob, Kon, rnd, forceAux;
	int i;
	double limitForce=0.06;


	for (i=0; i<CadhPack.ActiveCadh.size(); ++i )
	//for (int i=0; i<CadhPack.List.size(); ++i )
	{
		if (! CadhPack.List.at(i).full)
		{
		 	CadhPack.List.at(i).length= NodesDist( CadhPack.List.at(i).firstAdhPoint, 
		 												CadhPack.List.at(i).secondAdhPoint );
		 	CadhPack.List.at(i).force=CadhPack.GETEA()*(CadhPack.List.at(i).length- CadhPack.GETlengthBal())*CadhPack.List.at(i).cadherinNumber;
		 	forceAux=CadhPack.List.at(i).force/CadhPack.List.at(i).cadherinNumber;

		 	if (( forceAux < limitForce ) && (CadhPack.List.at(i).length > CadhPack.GETlengthBal()) )
		 	{
		 		Kon=CadhPack.GETKonRate()*CadhPack.GETdensity()*
		 			(1-(limitForce*150-forceAux)/(limitForce*150))*0.8;
		 		prob=(1-exp(-Kon*TimeStep))*0.75;
		 		rnd=distribution(generator);
		 		if (prob > rnd) 
		 		{
		 			CadhPack.Binds(i);
		 			Nodes.at(CadhPack.List.at(i).firstAdhPoint)-> SETbound(true);
					Nodes.at(CadhPack.List.at(i).secondAdhPoint)-> SETbound(true);
		 		}
	
			}
			else
			{}

		
		 	//(*it).force=(*it).length*CadhPack.GETK()/CadhPack.GETlengthBal();
	     }
	}

}

/*
void Monolayer::CheckReBinding()    //  Inverse Distance
{
	std::uniform_real_distribution<double> distribution(0, 1);
	double prob, Kon, rnd;
	int i;


	for (i=0; i<CadhPack.ActiveCadh.size(); ++i )
	//for (int i=0; i<CadhPack.List.size(); ++i )
	{
		if (! CadhPack.List.at(i).full)
		{
		 	CadhPack.List.at(i).length= NodesDist( CadhPack.List.at(i).firstAdhPoint, 
		 												CadhPack.List.at(i).secondAdhPoint );


		 	if ( CadhPack.List.at(i).length < CadhPack.GETbindLimit()*1.2 )
		 	{
		 		Kon=CadhPack.GETKonRate()*CadhPack.GETdensity()*
		 			(1-CadhPack.List.at(i).length/(CadhPack.GETbindLimit()*1.2));
		 		prob=(1-exp(-Kon*TimeStep));
		 		rnd=distribution(generator);
		 		if (prob > rnd) 
		 		{
		 			CadhPack.Binds(i);
		 			Nodes.at(CadhPack.List.at(i).firstAdhPoint)-> SETbound(true);
					Nodes.at(CadhPack.List.at(i).secondAdhPoint)-> SETbound(true);
		 		}
	
			}
			else
			{}

		}
		 	//(*it).force=(*it).length*CadhPack.GETK()/CadhPack.GETlengthBal();
	 }

}
*/
void Monolayer::CheckUnbinding()
{
	int ID, a, b, maxNum;
	double ubProbAux, rndNumber, forceAux,lengthAux;
	std::uniform_real_distribution<double> distribution(0, 1);

	for (auto it=CadhPack.ActiveCadh.begin(); it !=CadhPack.ActiveCadh.end(); )
	{
		ID=*it;
		a=CadhPack.List.at(ID).firstAdhPoint;
		b=CadhPack.List.at(ID).secondAdhPoint;
		lengthAux=NodesDist(a, b);
		CadhPack.List.at(ID).length=lengthAux;

//			forceAux=CadhPack.GETEA()/CadhPack.GETlengthBal()
//					*(CadhPack.List.at(ID).length- CadhPack.GETlengthBal());

		forceAux=CadhPack.GETEA()*(CadhPack.List.at(ID).length- CadhPack.GETlengthBal())*CadhPack.List.at(ID).cadherinNumber;								
		CadhPack.List.at(ID).force=forceAux;
		CadhPack.UpdateUbProb(ID);    // I divide for the number of cadherins in the cluster
		ubProbAux=CadhPack.List.at(ID).ubProb;
		maxNum=CadhPack.List.at(ID).cadherinNumber;
		for (int i=0; i < maxNum; i++)
		{
			rndNumber=distribution(generator);	
			if (rndNumber < ubProbAux) 
			{
				CadhPack.Unbinds(ID);
				if (CadhPack.List.at(ID).cadherinNumber == 0)
				{
					CadhPack.List.at(ID).active=false;
					it=CadhPack.ActiveCadh.erase(it);
					CadhPack.InactiveCadh.push_back(ID);
					UpdateMembpointFree(a,b,ID);
					CadhPack.List.at(ID).force=double(0);
					// Put memb points with bound fales and in the list
				}
				else if (i==maxNum-1)
					++it;
			}
			else if (i==maxNum-1)
				++it;
		}
	}

 	CadhPack.updateCts(TimeSlots);
}


void Monolayer::UpdateMembpointFree(int aIndex,int bIndex, int IDCadh)
{
	int aCellID, bCellID;



	if (Nodes.at(aIndex)-> GETFull())
	{
		aCellID=floor(aIndex/(nAdhPointCell+nCenterPointCell));  // Cell at which it is placed 
		Pool[aCellID].UnboundAdhPoint.insert(aIndex);
	}
	
	if (Nodes.at(bIndex)-> GETFull())
	{
		bCellID=floor(bIndex/(nAdhPointCell+nCenterPointCell));  // Cell at which it is placed 
		Pool[bCellID].UnboundAdhPoint.insert(bIndex);
	}

	
	Nodes.at(aIndex)-> SETCadhFree(IDCadh);
	Nodes.at(bIndex)-> SETCadhFree(IDCadh);

	//




//	std::cout<< " New Adh -> Cell: "<< aCellID << ";  Point: " << aIndex << std::endl;
//	std::cout<< "         -> Cell: "<< bCellID << ";  Point: " << bIndex << std::endl;


}

void Monolayer::CheckUnboundAdhVectors()   // when something binds again does not leave vector. This function fix this
{
	int ID;
	// Cell center
	for (int i=0; i < nCell; ++i)
	{
		for (auto it=Pool.at(i).UnboundAdhPoint.begin(); it != Pool.at(i).UnboundAdhPoint.end(); ) 
		{
			ID=*it;
			if (Nodes.at(ID)-> GETFull() == true)
			{
				it=Pool.at(i).UnboundAdhPoint.erase(it);
			}
			else
				it++ ;
		}
			 
	}

}


void Monolayer::InitialCheckofBoundSate()  // This is done to update availability of membrane points after initial conditions
{
	int aIndex;

	int aCellID;

	for (aIndex=0; aIndex < Nodes.size();aIndex++ ) 
	{
		if ((Nodes.at(aIndex)-> GETFull()) ==false)
		{
			aCellID=floor(aIndex/(nAdhPointCell+nCenterPointCell));  // Cell at which it is placed 
			Pool[aCellID].UnboundAdhPoint.insert(aIndex);
		}
	}	
}





void Monolayer::BalanceLangevin()
{

	double adaptiveTimeStep=TimeStep;
	int iter=1;
	double MaxXInc; 

	Eigen::VectorXd  FInt(2*meshSize);
//	Eigen::VectorXd  CInt(2*meshSize);
	Eigen::VectorXd  XFext(2*meshSize);
//	Eigen::VectorXd  Fvisc(2*meshSize);
//	Eigen::VectorXd  Xvisc(2*meshSize);
	Eigen::VectorXd  Faux(2*meshSize);
	Eigen::VectorXd  FRep(2*meshSize);

	FInt=Eigen::VectorXd::Zero(2*meshSize);
	XFext=Eigen::VectorXd::Zero(2*meshSize);
//	CInt=Eigen::VectorXd::Zero(2*meshSize);
//	Xvisc=Eigen::VectorXd::Zero(2*meshSize);
 	FRep=Eigen::VectorXd::Zero(2*meshSize);
 	Faux=Eigen::VectorXd::Zero(2*meshSize);
 	XIncLastStep=Eigen::VectorXd::Zero(2*meshSize);

 	std::vector<double> XIncIter;
 	XIncIter.resize(2*meshSize);

 	RepulsionMembraneGlobalInteraction(FInt);// Re Calculate closest points
	CalculateForceIntern(FInt);

//	RepulsionCellInCenter(FRep);

//	CalculateDamper(CInt);

//	Faux=FInt+ForceStep+ForceStepBr;
	Faux=FInt+ForceStep;
	for (int i=0; i < 2*meshSize; i++)
		XFext(i)=Faux(i)*TimeStep;
	//	XFext(i)=Faux(i)/(CIntGlob(i)+mediumViscCoef)*TimeStep;
	ApplyBdryConditions(XFext);
/*
	for (int i=0; i < 2*meshSize; i++)
		Fvisc(i)=(-1)*CInt(i)*XFext(i)/TimeStep;

	Xvisc=Fvisc/mediumViscCoef*TimeStep;
	ApplyBdryConditions(Xvisc);
*/	
	for (int i=0; i < 2*meshSize; i++)
		XIncIter[i]=XFext(i);
	
//	print_data_intern2(FInt,Faux , Fvisc,
//  						XFext, Xvisc, XIncIter, CInt);

	MaxXInc=*std::max_element(XIncIter.begin(),XIncIter.end());



	if ( MaxXInc > maximunDispAllowed )  // Make internal iterations
	{
		iter=floor(MaxXInc/maximunDispAllowed)+1;
		adaptiveTimeStep=TimeStep/iter;

		std::cout << "Number of iter:  " << iter << std::endl;
		/*
		for (int i=0; i < 2*meshSize; i++)
			XIncIter[i]=XIncIter[i]/iter;

		UpdateNodesPosition(XIncIter);
		for (int i=0; i < 2*meshSize; i++)
			XIncLastStep(i)+=XIncIter[i];
		*/

		for (int j=0; j<iter; ++j )
		{
			FInt=Eigen::VectorXd::Zero(2*meshSize);
			CalculateForceIntern(FInt);
	//		CalculateDamper(CInt);
			Faux=FInt+ForceStep+ForceStepBr;
			for (int i=0; i < 2*meshSize; i++)
				XFext(i)=Faux(i)*adaptiveTimeStep;
			//	XFext(i)=Faux(i)/(CIntGlob(i)+mediumViscCoef)*adaptiveTimeStep;
			ApplyBdryConditions(XFext);

//			for (int i=0; i < 2*meshSize; i++)
//				Fvisc(i)=(-1)*CInt(i)*XFext(i)/adaptiveTimeStep;
//
//			Xvisc=Fvisc/mediumViscCoef*adaptiveTimeStep;
//			ApplyBdryConditions(Xvisc);


	
			for (int i=0; i < 2*meshSize; i++)
				XIncIter[i]=XFext(i);

			UpdateNodesPosition(XIncIter);

			for (int i=0; i < 2*meshSize; i++)
			XIncLastStep(i)+=XIncIter[i];

//			print_data_intern2(FInt,Faux , Fvisc,
// 						XFext, Xvisc, XIncIter,XIncLastStep ,j, iter, CInt);

		}

		
	}

	else
	{
		for (int i=0; i < 2*meshSize; i++)
			XIncLastStep(i)+=XIncIter[i];

		UpdateNodesPosition(XIncLastStep);

		XIncTotal=XIncTotal+XIncLastStep;

	//	print_data_intern(FInt,Fvisc , XFext,
  	//					Xvisc);
	}

}

void Monolayer::ApplyBdryConditions(Eigen::VectorXd &VectorAux)
{
	int aux;
	for (int i=0; i < BdryIDs.size(); i++ )
	{
		aux=BdryIDs[i];
		VectorAux(aux)=double(0);
	} 

}
// To set bound bdry points
void Monolayer::SetBoundBdryConditions()
{
	int aux;
	for (int i=0; i < BdryIDs.size(); i++ )
	{
		aux=floor(BdryIDs[i]/2);
		Nodes.at(aux)-> SETBdry();
	} 

}


void Monolayer::CellRemodel()
{
	double distanceAux, cellAreaAux, remodelAux;
	int i,j, IDCenterAux, pointAux;
	
	for (j=0; j < Pool.size(); ++j)
	{

		// Remodel based on actual length

		cellAreaAux=double(0);
		IDCenterAux=Pool.at(j).GETIDCenter();
		for ( i=(IDCenterAux+1); i< (IDCenterAux+1+nAdhPointCell); ++i )
		{
			distanceAux=NodesDist(IDCenterAux, i);
			if (( distanceAux > aristLength*0.625)  && (distanceAux < aristLength*1.375 ))
			{
				remodelAux=Kremodel*(distanceAux-SFlength[i]);
				if (remodelAux >maximumRemodel)
					remodelAux=maximumRemodel;
				SFlength.at(i)=SFlength[i]+remodelAux;
				
			}
			cellAreaAux+=SFlength[i];
		}


		// Mantaining Area

		if (remodelUB)   // Only stress fibers that are unbound
		{

			remodelAux = double(-1)*(cellAreaAux-cellAreaInit)/double(Pool.at(j).UnboundAdhPoint.size()); // Romodel to each stress fiber unbound
			if (remodelAux > (maximumRemodel*remodelRatioVelVS) )
				remodelAux=maximumRemodel*remodelRatioVelVS;
			for (std::set<int>::iterator it=Pool.at(j).UnboundAdhPoint.begin(); it!=Pool.at(j).UnboundAdhPoint.end(); ++it)
			{
				pointAux=*it;
				distanceAux=NodesDist(IDCenterAux, pointAux);
				if (distanceAux < aristLength*1.375 )
					SFlength[pointAux]=SFlength[pointAux]+remodelAux;

	
			}
		}
		else   // Remodel to all stress fibers
		{
			remodelAux = double(-1)*(cellAreaAux-cellAreaInit)/double(nAdhPointCell); // Remodel to all stress fibers
			if (remodelAux >(maximumRemodel*remodelRatioVelVS))
				remodelAux=maximumRemodel*remodelRatioVelVS;
			for ( i=0; i< nAdhPointCell; ++i )
			{
				pointAux=IDCenterAux+1+i;
				distanceAux=NodesDist(IDCenterAux, pointAux);
				if (distanceAux < aristLength*1.5 )
					SFlength[pointAux]=SFlength[pointAux]+remodelAux;

			}
		}
/*

*/

	}
}

void Monolayer::Balance()
{
	std::ofstream myfile;
	Eigen::SimplicialCholesky<SpMat>  solver;  // performs a Cholesky factorization of A
 //	Eigen::ConjugateGradient<SpMat>  solver;  // performs a Cholesky factorization of A
 //	Eigen::SparseLU<SpMat>  solver;  // performs a Cholesky factorization of A

	Eigen::VectorXd  XIncIt(2*meshSize);
	Eigen::VectorXd  XIncStep(2*meshSize);
	Eigen::VectorXd  RNodes(2*meshSize);
	Eigen::VectorXd  FInt(2*meshSize);
	Eigen::VectorXd  FIntLast(2*meshSize);
	Eigen::VectorXd  FDamperIt(2*meshSize);
 	//	Eigen::VectorXd  FDamperStep(2*meshSize);

	Eigen::VectorXd  FDampertot(2*meshSize);

	SpMat KCadh(2*meshSize,2*meshSize);
	SpMat CCadh(2*meshSize,2*meshSize);
	SpMat KGlobal(2*meshSize,2*meshSize);
	SpMat CGlobal(2*meshSize,2*meshSize);

	FInt=Eigen::VectorXd::Zero(2*meshSize);
	XIncIt=Eigen::VectorXd::Zero(2*meshSize);
	RNodes=Eigen::VectorXd::Zero(2*meshSize);

	//	ForceGeneration();
	bool divergence=true;
	bool stop=false;
	int divider=1;
	int iter=0;
	int iterMax=8;
	double modRLast;
	double modR=1e10;
	double epsError=1e-4;

  while (divergence  && (!stop))
  {
  		iter=0;
	  	divergence=false;
	  	modR=1e10;

	  	XIncStep=Eigen::VectorXd::Zero(2*meshSize);	

	
		if (smallDefResolution==0) 
			BuildSparseMatrixMembBdry();
	
		BuilSparseMatrixCadhConnectivity(KCadh, CCadh);
	
		KGlobal= KCadh + KMemBdry;
		CGlobal= CCadh + CMemBdry;
	
		solver.compute(KGlobal+(CGlobal/TimeStep));
		if(solver.info()!=	Eigen::Success)
		 {
	  		std::cout<< "decomp failed" << std::endl;
	  		return;
	  	 }
	
	
		
		if (smallDefResolution==0)
		{ 
			FInt=Eigen::VectorXd::Zero(2*meshSize);
			CalculateForceIntern(FInt);
		}
		else
			CalculateForceInternAlternative(KGlobal,FInt,XIncTotal );
	//	RNodes= ForceInc + FInt-ForceInternLastStep;
	
		RNodes= ForceLastStep + FInt - FDamperAcum;
		
	
		ForceInc=Eigen::VectorXd::Zero(2*meshSize);   // **********************************************************************
		
	
		// Start iteration
		iter=0;
	 	modR=1e10;

	//    FDampertot=FDamperIt+FDamperAcum;

	// Uncoment to check convergence
	//    print_data_intern(KGlobal+(CGlobal/TimeStep),RNodes,FInt ,FDampertot,
	//  									XIncTotal, XIncIt, divider, iter);
		
	
	  	while( ( modR > epsError ) && (iter < iterMax) && (!divergence) ) 
	  	{
	  		
			
			// Calculate new increment
			RNodes=RNodes/divider;
			XIncIt=solver.solve(RNodes);
			if(solver.info()!=Eigen::Success) {
	  		std::cout<< "solving failed" << std::endl;
	  		return;
	  		}
			
			// Update position
			XIncStep=XIncStep+XIncIt;
			XIncTotal=XIncTotal+XIncIt;
			UpdateNodesPosition(XIncTotal);
	
			// Update for next iteration
	
			ForceInternLastStep=FInt;
			FInt=Eigen::VectorXd::Zero(2*meshSize);
	
			if (smallDefResolution==0) 
			{
				BuildSparseMatrixMembBdry();
	
				BuilSparseMatrixCadhConnectivity(KCadh, CCadh);
	
				KGlobal= KCadh + KMemBdry;
				CGlobal= CCadh + CMemBdry;
				solver.compute(KGlobal+(CGlobal/TimeStep));

				FInt=Eigen::VectorXd::Zero(2*meshSize);
				CalculateForceIntern(FInt);
			}
			else
				CalculateForceInternAlternative(KGlobal,FInt,XIncTotal );
	
			FDamperIt=(CGlobal/TimeStep)*(XIncStep);


			FDampertot=FDamperIt+FDamperAcum;
			RNodes= ForceLastStep + FInt - FDampertot;
	
	
	//			print_data_intern(KGlobal,RNodes,FInt , FDampertot,
	//  									XIncTotal, XIncIt, divider, iter);	
		/*	 
			 RNodes= ForceInc + FInt-ForceInternLastStep - ( (CGlobal/TimeStep)*(XIncIt-XIncLastStep) );
		*/	
	
	
			 XIncLastStep=XIncIt;
		
			 modRLast=modR;
			 modR=sqrt(RNodes.dot(RNodes));
		  	
		  	if (modR > modRLast) 
		  	{
		  	 	divergence=true; 	
	  			myfile.open ("error_convergence.txt");
	  			myfile << "DIVERGE step: " << step << ";   iter: " << iter << "\n";
	  			myfile.close();
	  			XIncTotal=XIncTotal-XIncStep;
	  			UpdateNodesPosition(XIncTotal);
	  			divider=divider*2;
	  			if (divider > 4000)
	  				stop=true;
	  			//iter=iterMax;
		  	 	//print_data_help(filopodium,pack, inc_x, ECM, loop_condition, it, step, convergence);
		  	 	//std::cout<< step << " :Divergencia !!!!!!!!!!!!!!!   Divergencia !!!!!!!!!!!!!!!" << std::endl;
		  	 }
	
		  	 
		  	 // print_residuo( R_nodes,inc_x,step,it);
		  	 ++iter;
	//	  	 std::cout<< "R mod: "<< sqrt(mod_R) << std::endl;
		  	 //std::cin.get();
		}
	
		
		// Guardar posicion inicial, funcion
	
		// Actualizar posiciones, funcion
    } // divergence

    if (modR > epsError)  
		{
			std::cout << step << ": Iteration max -*-*-*-*-*-*-*-*-*-*-*-*-: "<< modR << std::endl;
				std::cout << step << ": Divider: -*-*-*-*-*-*-*-*-*-*-*-*-:  "<< divider << std::endl;
			myfile.open ("error_convergence.txt");
	  		myfile << "step: " << step << ";   error: " << modR << "\n";
	  		myfile.close();
	  	}

	  	FDamperAcum=FDamperIt+FDamperAcum;
	
}
	

void Monolayer::UpdateCadhProp()
{
	int ID, CadhID;
	double EA;

	for (auto it=CadhPack.List.begin(); it !=CadhPack.List.end(); )
	{
		(*it).length=NodesDist( (*it).firstAdhPoint, (*it).secondAdhPoint );


	/*	 These properties have gone to the Adhesion Points
		(*it).gapCtAux=false;

		if ( (*it).active==false)	
			(*it).gapTotTime+=1.;
		else
			(*it).gapGenTime.at((*it).nGap)+=1;
	*/	
		it++;
	}

	for (int i=0; i<CadhPack.ActiveCadh.size();++i)
	{
		CadhID=CadhPack.ActiveCadh.at(i);
		EA=CadhPack.GETEA()*CadhPack.List.at(CadhID).cadherinNumber;
		CadhPack.List.at(CadhID).force= EA* (CadhPack.List.at(CadhID).length-CadhPack.GETlengthBal());
		
	}
}


void Monolayer::UpdateAdhPointCts()
{
	int ID;
 // Podria probvar este iterador. Chequear que va todo va en orden antes de hacerlo
 //	 for (auto& x: mymap) {
 //    std::cout << x.first << ": " << x.second << '\n';
 //  }

	for (auto i=0; i < Nodes.size(); i++)
	{
		
		Nodes.at(i)->gapCtAux=false;

		if ( Nodes.at(i)->get_bound()==false)	
			Nodes.at(i)->gapTotTime+=1.;
		else
			Nodes.at(i)->gapGenTime.at(Nodes.at(i)->nGap)+=1;
	}
}




void Monolayer::CalculateDamper(Eigen::VectorXd &CInt)
{
	int i, CadhID, centerID ;
	double viscCoef;
	for (i=0; i< ConnectivityMemb.at(0).size();++i)
		CalculateNodesDamper(CInt, membVisc, 
							ConnectivityMemb.at(0).at(i), ConnectivityMemb.at(1).at(i), false);
	

	for (i=0; i< ConnectivityBdry.at(0).size();++i)
		CalculateNodesDamper(CInt, bdryVisc, 
						ConnectivityBdry.at(0).at(i), ConnectivityBdry.at(1).at(i), true);



	// Add extra viscosity in nucleus
		for (i=0; i< Pool.size();++i)
		{
			centerID=Pool.at(i).GETIDCenter();
			CInt(2*centerID)+=fabs(nucleusViscCoef);	
			CInt(2*centerID+1)+=fabs(nucleusViscCoef);
		}



/*	for (i=0; i<CadhPack.ActiveCadh.size();++i)
	{
		CadhID=CadhPack.ActiveCadh.at(i);
		//viscCoef=CadhPack.GETvisc()*CadhPack.List.at(CadhID).cadherinNumber;
		viscCoef=CadhPack.GETvisc();		
		CalculateNodesDamper(CInt, viscCoef, CadhPack.List.at(CadhID).firstAdhPoint, 
								CadhPack.List.at(CadhID).secondAdhPoint, false);
	}
*/
}


void Monolayer::CalculateNodesDamper(Eigen::VectorXd &CInt, double viscCoef , int a, int b, bool bdry)
{
	double forcemod, lCurrent, EACurrent,lBalance;
	double sign;
	std::vector<double> VDir(2);

	lBalance=NodesDistInit(a, b);
	NodesVDir(a,b, VDir);

	lCurrent=NodesDist(a, b);
	forcemod=( lCurrent-lBalance);

	if (forcemod > 0)
		sign=double(+1);
	else
		sign=double(-1);

/*	
	CInt(2*b)+=viscCoef*VDir.at(0)*double(-1)*sign;
	CInt(2*b+1)+=viscCoef*VDir.at(1)*double(-1)*sign;
	if(!bdry)
	{
		CInt(2*a)+=viscCoef*VDir.at(0)*sign;
		CInt(2*a+1)+=viscCoef*VDir.at(1)*sign;
	}
*/

	CInt(2*b)+=fabs(viscCoef);	
	CInt(2*b+1)+=fabs(viscCoef);

	if(!bdry)
	{
		CInt(2*a)+=fabs(viscCoef);
		CInt(2*a+1)+=fabs(viscCoef);
	}

}


void Monolayer::CalculateForceIntern(Eigen::VectorXd &FInt)
{
	int i, CadhID ;
	double EA, KBendTotal, auxVisc;
	int cellIDAux;
	auxVisc=membVisc+ mediumViscCoef;
	for (i=0; i< ConnectivityMemb.at(0).size();++i)
	{	
		
		CalculateNodesForce(FInt, membEA, auxVisc, NodesDistInit( ConnectivityMemb.at(0).at(i),ConnectivityMemb.at(1).at(i) ), 
							ConnectivityMemb.at(0).at(i), ConnectivityMemb.at(1).at(i));
		CalculateBendingForce(FInt, KBendingMemb , mediumViscCoef, NodesAngleInit( ConnectivityMemb.at(0).at(i),ConnectivityMemb.at(1).at(i) ), 
							ConnectivityMemb.at(0).at(i), ConnectivityMemb.at(1).at(i));
	
	}
	auxVisc=bdryVisc+mediumViscCoef;
	for (i=0; i< ConnectivityBdry.at(0).size();++i)
	{
		
		CalculateNodesForceBdry(FInt, bdryEA, auxVisc, SFlength[ConnectivityBdry.at(1).at(i)] ,
						ConnectivityBdry.at(0).at(i), ConnectivityBdry.at(1).at(i));  
	// For SF repulsion
	//	cellIDAux=CellIDfromPoint(ConnectivityBdry.at(1).at(i));
	//	CalculateSFOverlapForce(FInt, ConnectivityBdry.at(0).at(i), ConnectivityBdry.at(1).at(i),cellIDAux);
	
		// Esta en false para que el centro tb se mueva 
	}

	auxVisc=mediumViscCoef;
	for (i=0; i<CadhPack.ActiveCadh.size();++i)
	{
		auxVisc=mediumViscCoef;
		CadhID=CadhPack.ActiveCadh.at(i);
		EA=CadhPack.GETEA()*CadhPack.List.at(CadhID).cadherinNumber;
		//CalculateNodesForce(FInt, EA, auxVisc, CadhPack.GETlengthBal(),
		//		CadhPack.List.at(CadhID).firstAdhPoint, CadhPack.List.at(CadhID).secondAdhPoint);
		// Repulsion
		CalculateNodesForceCadh(FInt, EA, auxVisc, CadhPack.GETlengthBal(),
				CadhPack.List.at(CadhID).firstAdhPoint, CadhPack.List.at(CadhID).secondAdhPoint);

		KBendTotal=CadhPack.GETKBending()*CadhPack.List.at(CadhID).cadherinNumber;
	//	CalculateBendingForce(FInt, KBendTotal, NodesAngleInit(CadhPack.List.at(CadhID).firstAdhPoint, CadhPack.List.at(CadhID).secondAdhPoint ) ,
	//			CadhPack.List.at(CadhID).firstAdhPoint, CadhPack.List.at(CadhID).secondAdhPoint);
	}

	//RepulsionCellInCenter(FInt);
	RepulsionMembraneGlobalForce(FInt);
	RepForceCenter(FInt);

}

void Monolayer::SFOverlapForceVTK(Eigen::VectorXd &FInt)
{
	int cellIDAux;
	for (int i=0; i< ConnectivityBdry.at(0).size();++i)
	{ 
		cellIDAux=CellIDfromPoint(ConnectivityBdry.at(1).at(i));
		CalculateSFOverlapForce(FInt, ConnectivityBdry.at(0).at(i), ConnectivityBdry.at(1).at(i),cellIDAux);
	}

}

void Monolayer::CalculateForceBendingMemb(Eigen::VectorXd &FInt)
{
	int i, CadhID ;
	double EA;
	for (i=0; i< ConnectivityMemb.at(0).size();++i)
		CalculateBendingForce(FInt, KBendingMemb , 1.,  NodesAngleInit( ConnectivityMemb.at(0).at(i),ConnectivityMemb.at(1).at(i) ), 
							ConnectivityMemb.at(0).at(i), ConnectivityMemb.at(1).at(i));

}

void Monolayer::CalculateForceBendingCadh(Eigen::VectorXd &FInt)
{
	int i, CadhID ;
	double  KBendTotal;
	
	for (i=0; i<CadhPack.ActiveCadh.size();++i)
	{
		CadhID=CadhPack.ActiveCadh.at(i);
		KBendTotal=CadhPack.GETKBending()*CadhPack.List.at(CadhID).cadherinNumber;
		CalculateBendingForce(FInt, KBendTotal, 1., NodesAngleInit(CadhPack.List.at(CadhID).firstAdhPoint, CadhPack.List.at(CadhID).secondAdhPoint ) ,
				CadhPack.List.at(CadhID).firstAdhPoint, CadhPack.List.at(CadhID).secondAdhPoint);
	}
}


void Monolayer::CalculateForceBdry(Eigen::VectorXd &FInt)
{
	int i, CadhID ;
	double EA;
	
	for (i=0; i< ConnectivityBdry.at(0).size();++i)
		CalculateNodesForce(FInt, bdryEA, 1., SFlength[ConnectivityBdry.at(1).at(i)],
						ConnectivityBdry.at(0).at(i), ConnectivityBdry.at(1).at(i));  
		// Esta en false para que el centro tb se mueva

}

void Monolayer::CalculateForceMemb(Eigen::VectorXd &FInt)
{
	int i, CadhID ;
	double EA;
	for (i=0; i< ConnectivityMemb.at(0).size();++i)
		CalculateNodesForce(FInt, membEA, 1.,NodesDistInit( ConnectivityMemb.at(0).at(i),ConnectivityMemb.at(1).at(i) ), 
							ConnectivityMemb.at(0).at(i), ConnectivityMemb.at(1).at(i));

}

void Monolayer::CalculateForceCadh(Eigen::VectorXd &FInt)
{
	int i, CadhID ;
	double EA;
	
	for (i=0; i<CadhPack.ActiveCadh.size();++i)
	{
		CadhID=CadhPack.ActiveCadh.at(i);
		EA=CadhPack.GETEA()*CadhPack.List.at(CadhID).cadherinNumber;
	//	CalculateNodesForce(FInt, EA, 1., CadhPack.GETlengthBal(),
	//						CadhPack.List.at(CadhID).firstAdhPoint, CadhPack.List.at(CadhID).secondAdhPoint);
		CalculateNodesForceCadh(FInt, EA, 1., CadhPack.GETlengthBal(),
							CadhPack.List.at(CadhID).firstAdhPoint, CadhPack.List.at(CadhID).secondAdhPoint);
	}
}
void Monolayer::CalculateAngleInit(Eigen::VectorXd &FInt)
{
	int i, a,b ;
	double  angleBalance;
//	angleBalance=AngleBalanceCenter(a,b);
	for (i=0; i< ConnectivityMemb.at(0).size();++i)
	{
		a=ConnectivityMemb.at(0).at(i);
		b=ConnectivityMemb.at(1).at(i);
		angleBalance=AngleBalanceCenter(a,b);
		FInt(2*a)+=cos(angleBalance);
		FInt(2*a+1)+=sin(angleBalance);
	}
}

void Monolayer::CalculateBendingForce(Eigen::VectorXd &FInt, double KBending, double barVisc,  double angleBalance , int a, int b)
{
	double forcemod, distX, distY, angleVal, angleValAux, distNodes;
	std::vector<double> VDir(2);

	
	angleVal=NodesAngle(a,b);

//	if ( angleVal < 0 ) 
//		angleValAux=angleVal+2*PI;


	angleBalance=AngleBalanceCenter(a,b);   // Bending stiffness like a circle. Comment for hexagon

	angleValAux=(angleVal - angleBalance);

	if ( angleValAux > PI ) 
		angleValAux=angleValAux-2*PI;
	else if (  angleValAux < (-PI ) )
		angleValAux=angleValAux+2*PI;

	forcemod=KBending*angleValAux/barVisc;

	// if ( ( (angleVal > PI/2.) && (angleVal > PI) ) || (angleVal < 0) && (angleVal >  (-2*PI/3) ) )
	//	forcemod=forcemod 


	// Arrange to avoid sf overlap



	distNodes=NodesDist(a, b);
	AngleVDir(angleVal-PI/2.,VDir );
	FInt(2*b)+=forcemod*VDir.at(0)/distNodes;
	FInt(2*b+1)+=forcemod*VDir.at(1)/distNodes;

	AngleVDir(angleVal+PI/2.,VDir );
	FInt(2*a)+=forcemod*VDir.at(0)/distNodes;
	FInt(2*a+1)+=forcemod*VDir.at(1)/distNodes;
}

void Monolayer::CalculateSFOverlapForce(Eigen::VectorXd &FInt , int centerP, int membP, int cellIDAux) 
{
	double forcemod, distX, distY, angle1, angle2, angleDiff,distNodes, angleVal, angleLimitAux;
	int membP2;
	bool RepulsionCond=false;
	std::vector<double> VDir(2);

	
	angle1=NodesAngle(centerP,membP);   // membP- centerP

	membP2=NormalizedID(membP+1, GETVertexID(cellIDAux,0));	
	angle2=NodesAngle(centerP,membP2);   // memb- center

	angleLimitAux=-0.002;   // Acts when it is closer to 25nm.    Angle with 240 bars is 0.026 rad between them. Put in negative since it is clockwise

	angleDiff=angle2-angle1;

	if ((angle1 < 0. ) && (angle2 > 0. ))
	{
		if (fabs(angleDiff) > PI/3.)
			RepulsionCond=false;
		else
			RepulsionCond=true;
	}
	else if ((angleDiff < angleLimitAux) && ( fabs(angleDiff) < PI/3.) ) 
		RepulsionCond=false;
	else
		RepulsionCond=true;


	// REst of 2-1 needs to be clockwise therefore negative, if difference is a very high negative value means that we are changing cuadrant

//	if ( angleVal < 0 ) 
//		angleValAux=angleVal+2*PI;

	if (RepulsionCond)
	{
		forcemod=(maxRadForce+maxCtForce)/(mediumViscCoef+bdryVisc)*0.85;
		AngleVDir(angle1+PI/2.,VDir );
		FInt(2*membP)+=forcemod*VDir.at(0);
		FInt(2*membP+1)+=forcemod*VDir.at(1);

		AngleVDir(angle2-PI/2.,VDir );
		FInt(2*membP2)+=forcemod*VDir.at(0);
		FInt(2*membP2+1)+=forcemod*VDir.at(1);
	}
		
	
}

double Monolayer::AngleBalanceCenter(int a, int b)
{
	double angleBalance;
	int auxCellID, centerID;
	//std::vector<double> middlePoint(2);
	
	//CalculateMiddlePoint(a,b,middlePoint);
	auxCellID=floor(a/(nAdhPointCell+nCenterPointCell));  // Cell at which it is placed 
	centerID=Pool.at(auxCellID).GETIDCenter();
	angleBalance=NodesAngle(centerID, a)-PI/double(2);
	return angleBalance;

}

bool Monolayer::CheckOldGap(GapType &auxGap)  // Return true if is equal to old gap. 
{
	bool GapNoRep=true;
	int i;
	i=0;
	while(  (i < GAPStepClosed.size() ) && GapNoRep )
	{	
		
		if ( GAPStepClosed.at(i).vertex == auxGap.vertex )  // If they are different there is no need ct print will update returning false
		{
			for (int j=0; j < auxGap.IDList.size();++j)
			{
				for (int z=0; z < GAPStepClosed.at(i).IDList.size();++z)
				{
					if (auxGap.IDList[j]== GAPStepClosed.at(i).IDList[z])
					{
						GapNoRep=false;
						auxGap.ctGapsID=GAPStepClosed.at(i).ctGapsID; 	// Copy the ct print

						std::cout << "RecicleGap: Init  " << auxGap.InitPoint << std::endl;
						std::cout << "odlGap: ID  " << GAPStepClosed.at(i).ctGapsID << std::endl;

					}
						
				}
			}

		}

		i++;	
	}
	return !GapNoRep;
}


void Monolayer::CheckGapFormation(int step)
{

	int cellID=floor(nCellRow/2);
	int startID=GETVertexID(cellID,0);
	std::vector<int> IDListAux;
	int ct=0;
	int ctAux=0;
	int initialPoint;
	int closestVertAux;
	double distVertexAux;
	double lengthAc=double(0);
	double acArea=double(0);
	double maxHeight=double(0);
	double averHeight=double(0);
	GapType auxGap;
	bool checkdirection; // In case the firs point is free to check the counter colcwise direction
	checkdirection=false;
	int ID;
	bool auxBool;

	ctAux=nAdhPointCell;
	if ( (Nodes.at(startID)->get_bound()==false) && (Nodes.at(startID)->GETGapIncluded()==false) )
	{	
		checkdirection=true;

		for (int i=0; i < nAdhPointCell; ++i)	// Check the counter clock wise direction
		{	
			ID=startID+ nAdhPointCell -i-1;
			ctAux=nAdhPointCell-i-1;
	
			if ( (Nodes.at(ID)->get_bound()==false) && (Nodes.at(ID)->GETGapIncluded()==false) )
			{
				initialPoint=ID;
				//lengthAc=sepAdhPoint*(ct+1);
				acArea+=MembDistance.at(ctAux).distance*sepAdhPoint;
				averHeight+=MembDistance.at(ctAux).distance;
				if (MembDistance.at(ctAux).distance > maxHeight)
					maxHeight=MembDistance.at(ctAux).distance;

				ct++;
				IDListAux.push_back(ID);
			}
			else  // Once there is something unbound lives the loop
			{
				i=nAdhPointCell;
			}
		}
	}

	ctAux=0;
	for (int i=startID; i < (startID+ nAdhPointCell); ++i)
	{
		if ( (Nodes.at(i)->get_bound()==false) && (Nodes.at(i)->GETGapIncluded()==false) )
		{
			if (ct==0)
				initialPoint=i;

			lengthAc=sepAdhPoint*(ct+1);
			acArea+=MembDistance.at(ctAux).distance*sepAdhPoint;
			averHeight+=MembDistance.at(ctAux).distance;
			if (MembDistance.at(ctAux).distance > maxHeight)
				maxHeight=MembDistance.at(ctAux).distance;

			IDListAux.push_back(i);
			ct++;
		}
		else
		{
			if (IDListAux.size() >= 1)
			{
				
				auxGap.startStep=step;
				auxGap.stepsActive=1;
				auxGap.Area=acArea;
				auxGap.InitPoint=initialPoint;
				auxGap.LastPoint=i-1;    // Minus one is good, trust yourself 
				auxGap.HeightAver=averHeight/double(ct);
				auxGap.HeightMax=maxHeight;
				auxGap.Length=lengthAc;
				auxGap.nAdhPoint=ct;
				auxGap.IDList=IDListAux;
				auxGap.active=true;
				CalculateGapArea(auxGap);

				if (auxGap.Area > gapThreshold)
				{
	
					//distvertex
	
					//Block point availables for gaps
					for (int j=0; j < ct; j++)
					{
						ID=IDListAux.at(j);
						(Nodes.at(ID))-> SETGapIncluded(true);
					}
					/*
					distVertexAux=MinumunDistanceToVertex(IDListAux, initialPoint, i, cellID, closestVertAux);
					auxGap.closestVert=closestVertAux;
					auxGap.DistVertex=distVertexAux;
					if (distVertexAux < vertexLimitDistance)
						auxGap.vertex=true;
					else
						auxGap.vertex=false;
					*/


					auxBool= CheckOldGap(auxGap);
					auxGap.IDGapType=GAP.size();

					if ( auxBool == false)	
					{
						if (auxGap.vertex)
						{
							if (ctVertexPrint==0)
								ctVertexFirst=step;
							auxGap.ctGapsID=ctVertexPrint;
							ctVertexPrint++;
						}
						else
						{
							if (ctInterfacePrint==0)
								ctInterfaceFirst=step;
		
							auxGap.ctGapsID=ctInterfacePrint;
							ctInterfacePrint++;
						}
						
					}
					
	
					GAP.push_back(auxGap);
					ActiveGAP.push_back(GAP.size()-1);
				}


			}
			lengthAc=double(0);
			averHeight=double(0);
			maxHeight=double(0);
			ct=0;
			acArea=double(0);
			IDListAux.clear();
		}

		ctAux++;
	}
}
void Monolayer::UpdateGaps(int step)
{
	// Last function
	int cellID=floor(nCellRow/2);
	int startIDRef=GETVertexID(cellID,0);
	int startID;
	std::vector<int> IDListAux;
	int itAux=0;
	int ctAux=0;
	int ct=0;
	int secondID;
	int initialPoint;
	int lastPointAux;
	int closestVertAux;
	double distVertexAux;
	double lengthAc=double(0);
	double acArea=double(0);
	double maxHeight=double(0);
	double averHeight=double(0);
	GapType auxGap;
	bool checkdirection=false; // In case the firs point is free to check the counter colcwise direction
	bool checkregap=false;
	bool vertexAux;
	bool overlapGap;

	int ID, IDGap;

	int stepsActiveAux;
	bool reduction=false; // To check if the gap is closing

//	for (auto it=GAP.begin(); it!=GAP.end();)
	GAPStepClosed.clear();
	for (int c=0; c <  ActiveGAP.size(); c++)
	{
	  overlapGap=true;	

	  IDGap=ActiveGAP.at(c);
	  if (GAP.at(IDGap).active)  // Needed to provent problems when a gap overlaps (it does not leave the vector until UpdateGAPVector() is called at the end)
	  {	
		checkregap=false;
		reduction=false;
		lengthAc=double(0);
		averHeight=double(0);
		maxHeight=double(0);
		ct=0;
		acArea=double(0);
		IDListAux.clear();

		stepsActiveAux=GAP.at(IDGap).stepsActive-1; //Check prevoius step
		itAux=0;

		for (int i=0; i < GAP.at(IDGap).IDList.size(); i++) // Check if the gap is closing
		{
			ID=GAP.at(IDGap).IDList.at(i);
			if (Nodes.at(ID)-> get_bound())
			{	
				//Nodes[ID]-> SETGapIncluded(false);
				reduction=true;
			}
		}

		if (reduction)
		{
			GAPStepClosed.push_back(GAP.at(IDGap));
			for (int i=0; i < GAP.at(IDGap).IDList.size(); i++) // Activate available for gap again to check
			{
				ID=GAP.at(IDGap).IDList.at(i);
				Nodes.at(ID)-> SETGapIncluded(false);
				GAP.at(IDGap).active=false;
				GAP.at(IDGap).endTime=step;
			}
			// Change comments if we want to start by closest to the vertex
			secondID=GAP.at(IDGap).LastPoint;		
			startID=GAP.at(IDGap).InitPoint;
            /*
			if (GAP.at(IDGap).closestVert==0)  // Vertex closest to start
			{
				secondID=GAP.at(IDGap).LastPoint;		
				startID=GAP.at(IDGap).InitPoint;
			}
			else
			{
				startID=GAP.at(IDGap).LastPoint;		
				secondID=GAP.at(IDGap).InitPoint;
			}
			*/
			ct=0;
			
			if ( (Nodes.at(startID)->get_bound()==false) && (Nodes.at(startID)->GETGapIncluded()==false) )
			{	
				

				for (int i=0; i < nAdhPointCell; ++i)	  // Counterclockwise
				{	
				
					ID=NormalizedID(startID-i-1, startIDRef);
					ctAux=ID-startIDRef;
					if ( (Nodes.at(ID)->get_bound()==false) && (Nodes.at(ID)->GETGapIncluded()==false) )
					{
						initialPoint=ID;
						lengthAc=sepAdhPoint*(ct+1);
						acArea+=MembDistance.at(ctAux).distance*sepAdhPoint;
						averHeight+=MembDistance.at(ctAux).distance;
						IDListAux.push_back(ID);
						ct++;
					}
					else  // Once there is something unbound lives the loop
					{
						i=nAdhPointCell;
					}

				}
				// Poner normal direction
			}

			// If Nodes beyond the last point are available for GAP in order to check later and go further in the loop if needed
			if ( (Nodes.at(secondID)->get_bound()==false) && (Nodes.at(secondID)->GETGapIncluded()==false) )   
			{
				for (int i=1; i < nAdhPointCell; ++i)	  // Clockwise
				{	
					ID=NormalizedID(secondID+i, startIDRef);
					ctAux=ID-startIDRef;
					if ( (Nodes.at(ID)->get_bound()==false) && (Nodes.at(ID)->GETGapIncluded()==false) )
					{
						itAux++;
					}
					else
						i=nAdhPointCell;

				}
			}

			for (int i=0; i < GAP.at(IDGap).IDList.size()+itAux; ++i)
			{
				ID=NormalizedID(startID+i, startIDRef);
				ctAux=ID-startIDRef;
				if ( (Nodes.at(ID)->get_bound()==false) && (Nodes.at(ID)->GETGapIncluded()==false) )
				{
					if (ct==0)
						initialPoint=ID;
	
					lengthAc=sepAdhPoint*(ct+1);
					acArea+=MembDistance.at(ctAux).distance*sepAdhPoint;
					averHeight+=MembDistance.at(ctAux).distance;
					if (MembDistance.at(ctAux).distance > maxHeight)
						maxHeight=MembDistance.at(ctAux).distance;
	
					IDListAux.push_back(ID);
					ct++;
				}
				else
				{
					if (acArea > gapThreshold*percentGapDissap)
					{
						checkregap=true;
						GAPStepClosed.pop_back();


						GAP.at(IDGap).stepsActive++;
						GAP.at(IDGap).Area=acArea;
						GAP.at(IDGap).InitPoint=initialPoint;
						GAP.at(IDGap).LastPoint=IDListAux.at(IDListAux.size()-1);
						GAP.at(IDGap).HeightAver=averHeight/double(ct);
						GAP.at(IDGap).HeightMax=maxHeight;
						GAP.at(IDGap).Length=lengthAc;
						GAP.at(IDGap).nAdhPoint=ct;
						GAP.at(IDGap).IDList=IDListAux;
						GAP.at(IDGap).active=true;

						//Check overlap

						overlapGap=CheckGapOverlap(GAP.at(IDGap) );
						while (overlapGap)
						{
							overlapGap=CheckGapOverlap(GAP.at(IDGap) );
						}


						//distvertex
						if (IDListAux.size()!= ct)
							std::cout<< "Error size check 1768 line" << std::endl;


						//Block point availables for gaps
						for (int j=0; j < IDListAux.size(); j++)
						{
							ID=IDListAux.at(j);
							(Nodes.at(ID))-> SETGapIncluded(true);
						}
						distVertexAux=MinumunDistanceToVertex(IDListAux, initialPoint, i, cellID, closestVertAux);
						GAP.at(IDGap).closestVert=closestVertAux;
						GAP.at(IDGap).DistVertex=distVertexAux;
						
						/*
						if (distVertexAux < vertexLimitDistance)
							auxGap.vertex=true;
						else
							auxGap.vertex=false;
						*/
						vertexAux=GAP.at(IDGap).vertex;
						CalculateGapArea(GAP.at(IDGap));
						if (GAP.at(IDGap).vertex != vertexAux)
						{
							if (GAP.at(IDGap).vertex)
							{
								GAP.at(IDGap).ctGapsID=ctVertexPrint;
								ctVertexPrint++;
							}
							else
							{
								GAP.at(IDGap).ctGapsID=ctInterfacePrint;
								ctInterfacePrint++;
							}
						}



						i = GAP.at(IDGap).IDList.size()+itAux; //End the loop
					}
	
					lengthAc=double(0);
					averHeight=double(0);
					maxHeight=double(0);
					ct=0;
					acArea=double(0);
					IDListAux.clear();
					itAux=0;
				}
	
				
			}
			// Check if the there is still gap big enough to be, starting by the closest part to vertex
			// Break it
		}
		else // Check if the gap is increasing
		{
		//	overlapGap=CheckGapOverlap(GAP.at(IDGap) );// 
			checkregap=true;
		    overlapGap=true;	
		    //	std::cout<< "increase" << std::endl;
		  	while (overlapGap)
		  	{

				IDListAux=GAP.at(IDGap).IDList; // Copy from the previous step
				lastPointAux=GAP.at(IDGap).LastPoint;
				initialPoint= GAP.at(IDGap).InitPoint;
	
				startID=GAP.at(IDGap).InitPoint;
				for (int i=1; i < nAdhPointCell; ++i)	  // Counterclockwise
				{	
						
					ID=NormalizedID(startID-i, startIDRef);
			
					if ( (Nodes.at(ID)->get_bound()==false) && (Nodes.at(ID)->GETGapIncluded()==false) )
					{
						initialPoint=ID;
						IDListAux.push_back(ID);
					}
					else  // Once there is something unbound lives the loop
					{
						i=nAdhPointCell;
					}
				}
	
				secondID=GAP.at(IDGap).LastPoint;
				for (int i=1; i < nAdhPointCell; ++i) //clockwise
				{
					ID=NormalizedID(secondID+i, startIDRef);
					if ( (Nodes.at(ID)->get_bound()==false) && (Nodes.at(ID)->GETGapIncluded()==false) )
					{
						IDListAux.push_back(ID);
						lastPointAux=ID;
					}
					else  // Once there is something unbound lives the loop
					{
						i=nAdhPointCell;
					}
				}
						
				// Update GAP
				 GAP.at(IDGap).IDList=IDListAux;
				 GAP.at(IDGap).nAdhPoint=IDListAux.size();
				 GAP.at(IDGap).InitPoint=initialPoint;
				 GAP.at(IDGap).LastPoint=lastPointAux;
				 GAP.at(IDGap).Length=IDListAux.size()*sepAdhPoint;
				 distVertexAux=MinumunDistanceToVertex(IDListAux, initialPoint, lastPointAux, cellID, closestVertAux);
				 GAP.at(IDGap).closestVert=closestVertAux;
				 GAP.at(IDGap).DistVertex=distVertexAux;
							
				// Check overlap!
				overlapGap=CheckGapOverlap(GAP.at(IDGap) );
			
		    } // End overlap


			for (int i=0; i< IDListAux.size(); i++)
			{
				ID=IDListAux[i];
				ctAux=ID-startIDRef;

				(Nodes.at(ID))-> SETGapIncluded(true);  // Set nodes unavialable for new gap

				acArea+=MembDistance.at(ctAux).distance*sepAdhPoint;
				averHeight+=MembDistance.at(ctAux).distance;
				if (MembDistance.at(ctAux).distance > maxHeight)
					maxHeight=MembDistance.at(ctAux).distance;

			}

			GAP.at(IDGap).Area=acArea;
			GAP.at(IDGap).HeightAver=averHeight/double(IDListAux.size());
			GAP.at(IDGap).HeightMax=maxHeight;
			//Remove if we want to change the label vertex or not of the gap during time
			/*  
			if (distVertexAux < vertexLimitDistance)
				GAP.at(IDGap).vertex=true;
			else
				GAP.at(IDGap).vertex=false;
			*/
			GAP.at(IDGap).stepsActive++;
			

			vertexAux=GAP.at(IDGap).vertex;
			CalculateGapArea(GAP.at(IDGap));
			if (GAP.at(IDGap).vertex!= vertexAux)  // if its condition of vertex/ interface changes update in order to print
			{
				if (GAP.at(IDGap).vertex)
				{
					GAP.at(IDGap).ctGapsID=ctVertexPrint;
					ctVertexPrint++;
				}
				else
				{
					GAP.at(IDGap).ctGapsID=ctInterfacePrint;
					ctInterfacePrint++;
				}
			}
		}
	/*	
		if	(checkregap)
		{
			it++;
		}
		else
		{
			GAP.at(IDGap).active=false;
			GAP.at(IDGap).endTime=step;
			InactiveGAP.push_back(GAP.at(IDGap));
			it=GAP.erase(it);
		}
	*/



	  } // En if	
	}
	UpdateGAPVector();
}

void Monolayer::UpdateGAPVector()
{
	ActiveGAP.clear();
	InactiveGAP.clear();
	for (int i=0; i< GAP.size(); ++i)
	{
		if (GAP.at(i).active)
			ActiveGAP.push_back(i);
		else
			InactiveGAP.push_back(i);
	}
}

void Monolayer::CalculateGapArea(GapType &GapAux)
{

	int ctCell=1;
	int ct=0;
	int  IDAux, IDChosen,cellIDChosen;
	int ctJump=0;
	int ctAux=0;;
	std::vector<int> boundListAux;
	int pointPrevCell;
	int	IDJumpChosen, IDAuxFirst;
	int IDJump;

	int maxAux;
	int minAux;
	//std::vector<int> cellIDList;
	GapAux.CurrentShapeIDList.clear();
	
	int startP=GapAux.InitPoint;
	int cellStartP, startAux;
	int ID=startP;
	int cellIDAux=CellIDfromPoint(startP);
	//cellIDList.push_back(cellIDAux);
	startAux=startP;
	int boundListAuxSize;

	int lastCellID, lastPointID;  // Cell and point that binds to the previous one to initial point 
	int prevCellID, prevPointID;
	
	bool backToWhile;
	double angle, angleAux;
	std::vector<double> vDirAux1, vDirAux2;
	vDirAux1.resize(2);
	vDirAux2.resize(2);

	std::vector<std::vector<int>> RepShapeID;
	std::vector<bool> RepVertex;
	int repetition;
	bool secondLevel, exitLevel, exitThirdLevel;	
	std::vector<int> ListCellNotGoBack;
	int ctLevel=0;



	IDAux=NormalizedID(startP-1, GETVertexID(cellIDAux,0));   // Chack previous point to the start to check that it is bound
	
	  // Vector Size it is always maximun, if there is a -1 means that this side is not bound
	  //	if (minAux < 0 )
	  //		boundListAuxSize=1;
	  //	else
	  //		boundListAuxSize=2;

	GetBoundListOK(boundListAux, IDAux); 

	repetition=boundListAux.size();
	RepShapeID.resize(repetition);
	RepVertex.resize(repetition);
	repetition=1;

  	for  (int rep=0; rep < repetition; rep++)
 	{	
  		ID=startP;
		cellIDAux=CellIDfromPoint(startP);
  		IDAuxFirst=NormalizedID(startP-1, GETVertexID(cellIDAux,0));   // Check previous point to the start to check that it is bound
  		GetBoundListOK(boundListAux, IDAux); 

  		lastPointID=boundListAux[rep];
		lastCellID=CellIDfromPoint(lastPointID);  // To resolve doublts when cell has two options for moving 
		// 1st Level of decission: Find the first point bound and take it. Try to find previous cell, otherwise go to the other cell
 		ctCell=1;

 		ListCellNotGoBack.clear();
 		//ct=0;
 		secondLevel=true;
		// Poner IDAux ¿?
 		while (secondLevel)
 		{
 			ctLevel++;
 		//	std::cout << "Second level: " <<ctLevel		 <<std::endl;

			exitLevel=false;
			ct=0;
 			while(!exitLevel)   // Dont exit unless we find a jum that leads us to a cell we haven been there before
 			{	
 					
				while (Nodes.at(ID)->get_bound() == false )
				{
					RepShapeID.at(rep).push_back(ID);
					cellIDAux=CellIDfromPoint(ID);
					ID++;
					ID=NormalizedID(ID, GETVertexID(cellIDAux,0));
					
					if (ct > 3*nAdhPointCell)
						abort_program(); 
					ct++;
				} 
				
			//	JumpCell(ID, cellIDAux, lastCellID);
				GetBoundListOK(boundListAux, ID);
				IDJumpChosen=boundListAux[0];
				cellIDAux=CellIDfromPoint(IDJumpChosen);

				exitLevel=true;
				for(int i=0; i< ListCellNotGoBack.size(); i++) //Check if we go back to cell we have been before. Dont do that
				{
					if (ListCellNotGoBack.at(i)==cellIDAux)
						exitLevel=false;
				}
				if (exitLevel==false)
				{
					RepShapeID.at(rep).push_back(ID);
					cellIDAux=CellIDfromPoint(ID);
					ID++;
					ID=NormalizedID(ID, GETVertexID(cellIDAux,0));
				//	std::cout << "REturn to level" << std::endl;
				}
				else  // If we jum to other cell check that if we go back to same cell we are we are not going to hingher position
				{
					if (cellIDAux==prevCellID)  // We go back to cell we where (only posible in first while iteration)
					{
						if (IsMajorInCell(IDJumpChosen,prevPointID))  // True if first > second
						{
							RepShapeID.at(rep).push_back(ID);
							cellIDAux=CellIDfromPoint(ID);
							ID++;
							ID=NormalizedID(ID, GETVertexID(cellIDAux,0));
						//	std::cout << "REturn to level" << std::endl;
							exitLevel=false;		

						}

					}
				}
			}
			

			IDAux=ID;
			IDChosen=ID;

			cellIDAux=CellIDfromPoint(ID);
			GetBoundListOK(boundListAux, ID);
// I use this for the starting case where there is point bound to two cells
			if (boundListAuxSize>1)
			{

				if (CellIDfromPoint(boundListAux[0]) == CellIDfromPoint(boundListAux[1])) // If both Point are form the same cell choose the biggest
				{
					IDJumpChosen=ReturnMajorInCell(boundListAux[0], boundListAux[1]);
				}
				else	//Go to lastcellID
				{
					if( CellIDfromPoint(boundListAux[0])==lastCellID ) 
						IDJumpChosen=boundListAux[0];
					else if	( CellIDfromPoint(boundListAux[1])==lastCellID ) 			
						IDJumpChosen=boundListAux[1];	
					else // If non of the are lastcellID
					{
						IDJumpChosen=boundListAux[0];	// I put it so it continuous but it is wrong check if it happens
					//	std::cout << "ERROR: Check first level final. GAP LOST" << std::endl;
					}
				}	
			
	 		}
	 		else
	 			IDJumpChosen=boundListAux[0];



		


		/*	for (int ii=1; ii<2; ++ii)
			{
				cellIDAux=CellIDfromPoint(ID);
				ID=NormalizedID(IDAux+ii, GETVertexID(cellIDAux,0));
				if ( Nodes.at(ID)->get_bound() == true )
				{
					GetBoundListOK(boundListAux, ID);
					IDJump=boundListAux[0];
					if ( CellIDfromPoint(IDJump) ==	CellIDfromPoint(IDJumpChosen) )  // Same cell, We chose the maximun numbler
	 				{
	 					if (IDJump > IDJumpChosen)
	 					{
	 						IDChosen=ID;
	 						IDJumpChosen=IDJump;
	 					}
	
	 				}
	 				else  // Diferent cell make the angle thing or leave it like it is
	 				{
	
	 				}
	
				 }
			 }
			
			*/
	
			// PUT ID is bound 
			
			startAux=IDJumpChosen; // In order to check angles
			prevPointID=IDChosen;
			prevCellID=CellIDfromPoint(IDChosen);
	
			RepShapeID.at(rep).push_back(IDChosen);
			RepShapeID.at(rep).push_back(IDJumpChosen);		//Point of new cell added
			ID=IDJumpChosen;
			ctCell++;
			cellIDAux=CellIDfromPoint(ID);
			ListCellNotGoBack.push_back(cellIDAux);

			// Check if we are in LasCellID
		 	if ( cellIDAux == lastCellID )
		 	{ 
		 		secondLevel=false;
		 	}
		 	else
		 	{	
		 		secondLevel=true;
		 		ID++;
				ID=NormalizedID(ID, GETVertexID(cellIDAux,0));
			}
		}// END SECOND LEVEL
		// Try to find lastcellID, otherwise go to the other cell
/*		if (boundListAuxSize>1)
		{

			if (CellIDfromPoint(boundListAux[0]) == CellIDfromPoint(boundListAux[1])) // If both Point are form the same cell choose the biggest
			{
				ID=ReturnMajorInCell(boundListAux[0], boundListAux[1]);
			}
			else	//Go to lastcellID
			{
				if( CellIDfromPoint(boundListAux[0])==lastCellID ) 
					ID=boundListAux[0];
				else if	( CellIDfromPoint(boundListAux[1])==lastCellID ) 			
					ID=boundListAux[1];	
				else // If non of the are lastcellID
				{
					ID=boundListAux[1];	// I put it so it continuous but it is wrong check if it happens
					std::cout << "ERROR: Check first level final. GAP LOST" << std::endl;
				}
			}	
		
	 	}
	 	else
	 		ID=boundListAux[0];

*/
	 


	 	// Go to second level
		// In second level we go to all the points until we reach a point that it is bound to LastCallID


	 	// We are in third level now:
	 	// Check that we did not pass the lastPointID
	 	// Go to the nodes until we are in lastPointID
	 	ct=0;
	 	if (IsMajorInCell(ID,lastPointID))
	 		exitThirdLevel=true;  // Point is already on the List
		else
		{
			ID++;
			ID=NormalizedID(ID, GETVertexID(cellIDAux,0));
			exitThirdLevel=false;
		}

		while (!exitThirdLevel)
		{
		//	std::cout << "Exit exitThirdLevel"<< std::endl;
			RepShapeID.at(rep).push_back(ID);
		
			if (ID==lastPointID)
				exitThirdLevel=true;
			else
			{
				ID++;
				ID=NormalizedID(ID, GETVertexID(cellIDAux,0));
				
				if (ct > 3*nAdhPointCell)
					abort_program();
				ct++;
			}		
		}	

		RepShapeID.at(rep).push_back(IDAuxFirst);  // Add point previous to initial one (startP-1)

		// Check if it is vertex
		if (ctCell<=2)
			RepVertex.at(rep)=false;
		else
			RepVertex.at(rep)=true;


	}

	double area0, area1;
	area0=AreaPoligonSymple(RepShapeID.at(0));
	area0=polygonArea(RepShapeID.at(0));

	if (repetition > 1) // We chose the minimun area
	{
		
	//	area1=AreaPoligonSymple(RepShapeID.at(1));	
		area1=polygonArea(RepShapeID.at(1));

		if ( area0 < area1)
		{
			GapAux.Area=area0;
			GapAux.vertex=RepVertex.at(0);
			GapAux.CurrentShapeIDList=RepShapeID.at(0);
		}
		else
		{
			GapAux.Area=area1;
			GapAux.vertex=RepVertex.at(1);
			GapAux.CurrentShapeIDList=RepShapeID.at(1);	
		}

	}
	else
	{
		GapAux.Area=area0;
		GapAux.vertex=RepVertex.at(0);
		GapAux.CurrentShapeIDList=RepShapeID.at(0);
	}

	// Temporary To avoid crazy cells
	//if (GapAux.Area > cellRealArea/2.2)
	//	GapAux.Area=0.;

	// Avoid crazy cells
	int ctBound;
	ctBound=0;
	for (int jj=0;jj< GapAux.CurrentShapeIDList.size(); jj++)
	{
		ID=GapAux.CurrentShapeIDList.at(jj);
		if (Nodes.at(ID)->get_bound() == true)
			ctBound++;

	}

	if ( ctBound > floor(nAdhPointSide/2.6) )
		GapAux.Area=0.;

} 



void Monolayer::abort_program()
{
	char name[100];
	char orden[1024]; 
	
	FILE*fOut;
	sprintf(name, "FAILURE_GAP_%s.txt" , index_cluster_Aux);			// Number of Gaps generated
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	fclose(fOut);
	exit (EXIT_FAILURE);
}
/*

void Monolayer::CalculateGapArea(GapType &GapAux)
{
	int ctCell=1;
	int ct=0;
	int  IDAux;
	int ctJump=0;
	int ctAux=0;;
	std::vector<int> boundListAux;
	int pointPrevCell;

	int maxAux;
	int minAux;
	std::vector<int> cellIDList;
	GapAux.CurrentShapeIDList.clear();
	
	int startP=GapAux.InitPoint;
	int cellStartP;
	int ID=startP;
	int cellIDAux=floor(startP/(nAdhPointCell+1));
	cellIDList.push_back(cellIDAux);

	int boundListAuxSize;

	int previousCellID=cellIDAux; // Previous cell we have gone trugh
	int lastCellID, lastPoint;  // Cell and point that binds to the previous one to initial point 
	
	bool loopClosed=false;
	bool backToWhile;
	bool quickJump=false;
	bool sameCellCondition,sameCellConditionBacwards;


	std::vector<std::vector<int>> RepShapeID;
	std::vector<bool> RepVertex;
	int repetition;



	IDAux=NormalizedID(startP-1, GETVertexID(cellIDAux,0));   // Chack previous point to the start to check that it is bound
	//boundListAux=Nodes.at(IDAux)-> CadhPointer;
	GetBoundListOK(boundListAux, IDAux); 

	repetition=boundListAux.size();
	RepShapeID.resize(repetition);
	RepVertex.resize(repetition);


  for  (int rep=0; rep < repetition; rep++)
  {	
  	lastPoint=boundListAux[i];
	lastCellID=floor(boundListAux[i]/(nAdhPointCell+1));  // To resolve doublts when cell has two options for moving 
	// 1st Level of decission: Find the first point bound and take it. Try to find previous cell, otherwise go to the other cell
 
	while (!loopClosed)
	{
		//std::cout << "LoopClosed enter  "  << std::endl;
		ct=0;
		backToWhile=true;
		//crazyCondition= ( (ct > 0) && () )
	  while (backToWhile)
	  {
	  	 backToWhile=false;
		while (Nodes.at(ID)->get_bound()== false )
		{
			RepShapeID.at(rep).push_back(ID);
			cellIDAux=floor(ID/(nAdhPointCell+1));
			ID++;
			ID=NormalizedID(ID, GETVertexID(cellIDAux,0));
			ct++;
		} 

		
		if (ct < 0)  // Si al siguiente de cambiarte quieres volver a la mima célula lo evita   Antes estaba en 2, chequear si el cambio a 1 usimepre esta bien
		{
			GetBoundListOK(boundListAux, ID);
			boundListAuxSize=boundListAux.size();
			ctAux=0;
			for (int i=0; i<boundListAuxSize; i++)
			{

				sameCellConditionBacwards=( ( ((-pointPrevCell+boundListAux[i])>0) && ( (-pointPrevCell+boundListAux[i]) < nAdhPointSide ) ) 
				  || ( ((-pointPrevCell+nAdhPointCell+boundListAux[i])>0) && ( (-pointPrevCell+nAdhPointCell+boundListAux[i]) < nAdhPointSide) ) );
	
				if ((floor(boundListAux[i]/(nAdhPointCell+1)) == previousCellID ) && sameCellConditionBacwards )  // No creo que haga falta chequer si hay dos, si los hay debería ir a la misma celula
				{	// Only jump when goes to the point of the same cell as prevoius one in a colser UpdateNodesPosition
					backToWhile=true;
					if (ctAux==0)
					{	
						RepShapeID.at(rep).push_back(ID);
						cellIDAux=floor(ID/(nAdhPointCell+1));
						ID++;
						ID=NormalizedID(ID, GETVertexID(cellIDAux,0));
						ctAux++;
					}	
				}
				else
				{
					backToWhile=false;
					i=boundListAuxSize-1;
				}
				
			}
		}
	  }
		
	  // PUT ID is bound 
	  RepShapeID.at(rep).push_back(ID);
	  //boundListAux=Nodes.at(ID)-> CadhPointer;   // Lo tenia comentado, no se poque
	  GetBoundListOK(boundListAux, ID);
	  boundListAuxSize=boundListAux.size();
		
	  // Size it always maximun, if there is a -1 means that this side is not bound
	  //	if (minAux < 0 )
	  //		boundListAuxSize=1;
	  //	else
	  //		boundListAuxSize=2;

		// Change to new cell, choose between the options the one that correspondes to the maximun value

		// Hay que chequear que la no este unido a dos celulas diferentes. Es facil, mirar el punto anterior al primer punto a la celula a la que esta unido y coger la misma en caso de duda
		// preparar mensaje de error con el step en que se da
		// IF TWO DIFFERENT CELLS

		if (boundListAuxSize>1)
		{
			maxAux= *max_element(boundListAux.begin(), boundListAux.end());
			minAux= *min_element(boundListAux.begin(), boundListAux.end());
			// 
			if ((floor(boundListAux[0]/(nAdhPointCell+1)) != floor(boundListAux[1]/(nAdhPointCell+1)) ) )  // Check if it is bound to more than one cell
			{	
				if (ctCell==1)  // First time: look for contrary cell
				{
					if (floor(boundListAux[0]/(nAdhPointCell+1)) == lastCellID )  // Choose the cell closest, avoid vertex and become interface
						ID=boundListAux[0];
					else
						ID=boundListAux[1];
				}
				else// Second time look for original cell
				{
					if ( (floor(boundListAux[0]/(nAdhPointCell+1)) == cellIDList[0] ) && (floor(boundListAux[1]/(nAdhPointCell+1)) == cellIDList[0] ) )
					{
						if (fabs(maxAux-minAux) > floor(nAdhPointCell/2))
							ID = minAux;
						else
							ID = maxAux;
					}

					else if (floor(boundListAux[0]/(nAdhPointCell+1)) == cellIDList[0] )  // Choose the cell closest, avoid vertex and become interface
						ID=boundListAux[0];
					else
						ID=boundListAux[1];

				}
			}	
			else
			{
				if (fabs(maxAux-minAux) > floor(nAdhPointCell/2))
					ID = minAux;
				else
					ID = maxAux;
			}
			
		}
		else
			ID = boundListAux[0];	
		

		// ADD NEW POINT OF NEW CELL
		


		pointPrevCell=GapAux.CurrentShapeIDList[RepShapeID.at(rep).size()-1];
		RepShapeID.at(rep).push_back(ID);
		cellIDAux=floor(ID/(nAdhPointCell+1));
		ctCell++;
		previousCellID=cellIDList[cellIDList.size()-1];
		cellIDList.push_back(cellIDAux);
		GetBoundListOK(boundListAux, ID);
		boundListAuxSize=boundListAux.size();
		

		if (cellIDAux==cellIDList[0])
			loopClosed=true;
		else
		{

			// IF IT IS BOUND TO A DIFFERENT CELL THAN THE ONE THAT WE ARE COMING FROM, GO THERE
			if (boundListAuxSize>1)
			{
	
				for (int i=0;i < boundListAux.size();i++)
				{
					sameCellCondition=( ( ((pointPrevCell-boundListAux[i])>0) && ( (pointPrevCell-boundListAux[i]) < nAdhPointSide ) ) 
					|| ( ((pointPrevCell+nAdhPointCell-boundListAux[i])>0) && ( (pointPrevCell+nAdhPointCell-boundListAux[i]) < nAdhPointSide) ) );
					// To ensure that if it goes there, the point is in the range of the cell where we come form
					if (( floor(boundListAux[i]/(nAdhPointCell+1)) != cellIDList[cellIDList.size()-2] )|| sameCellCondition) // Check if it is bound to different cell that the one we come from
					{	
						quickJump=true;
						ID=boundListAux[i]; // Change to this new cell
					}	
				}
			}
	
			if (quickJump) // CANNOT BE MORE THAN A QUICKJUMP IN A ROW, IT COULD NOT BE GAP IF THIS HAPPENS
			{
				RepShapeID.at(rep).push_back(ID);
				cellIDAux=floor(ID/(nAdhPointCell+1));
				ctCell++;
				cellIDList.push_back(cellIDAux);
				quickJump=false;
				if (cellIDAux==cellIDList[0])
					loopClosed=true;
			}
			else  
			{
				ID++;
				ID=NormalizedID(ID, GETVertexID(cellIDAux,0));
			}
			
		}
		if (ctCell >4)
		{
			loopClosed=true;
			std::cout << "LoopClosed exited because cells involded > 3  "  << std::endl;
			generalRupture=true;
		}
		
	}

	ctCell--; // Remove the repeated one
	// Real Vertex
	if (ctCell<=2)
		GapAux.vertex=false;
	else
		GapAux.vertex=true;



   }
	//GapAux.Area=polygonArea(GapAux.CurrentShapeIDList);
	GapAux.Area=AreaPoligonSymple(GapAux.CurrentShapeIDList);
	
	// aNALYZE IF VERTEX

} 
*/

bool Monolayer::IsPoligonSymple(std::vector<int> PointList )
{
	std::vector<PointCGal> PolPoints;
	PolPoints.reserve(PointList.size());
  	int ID;
  	bool simplePol;

  	for (int i=0; i<PointList.size(); i++) 
  	{
  		ID=PointList[i];
  		PolPoints.push_back(PointCGal(Nodes.at(ID)->Coord[0],Nodes.at(ID)->Coord[1]) );
    }

    Polygon_2 pgn(PolPoints.begin(), PolPoints.end());

    simplePol=pgn.is_simple();
  	return simplePol; 
}

double Monolayer::AreaPoligonSymple(std::vector<int> PointList )
{
	std::vector<PointCGal> PolPoints;
	PolPoints.reserve(PointList.size());
  	int ID;
  	bool AreaPol;

  	for (int i=0; i<PointList.size(); i++) 
  	{
  		ID=PointList[i];
  		PolPoints.push_back(PointCGal(Nodes.at(ID)->Coord[0],Nodes.at(ID)->Coord[1]) );
    }

    Polygon_2 pgn(PolPoints.begin(), PolPoints.end());

    AreaPol=pgn.area();
  	return AreaPol; 

}


// Index rep is the position in activeGAP of &GapAux
bool Monolayer::CheckGapOverlap(GapType &GapAux)
{
	int cellID=floor(nCellRow/2);
	int startIDRef=GETVertexID(cellID,0);	
	int IDGap, ID;
	int initRef=GapAux.InitPoint;
	int lastRef=GapAux.LastPoint;
	bool overlapGap=false;
	for (int c=0; c <  ActiveGAP.size(); c++)
	{
	  IDGap=ActiveGAP.at(c);

// 	  if (GAP.at(IDGap).active)  // When I do it recursive in case of three gaps I have to add this to not get the same one all the time	
// 	  { 	
		if ((GAP.at(IDGap).IDGapType != GapAux.IDGapType) && GAP.at(IDGap).active)// it is the same 
		{
			
			ID=NormalizedID(GAP.at(IDGap).InitPoint-1, startIDRef);
			if (lastRef == ID)  // Se unen
			{
				GapAux.LastPoint=GAP.at(IDGap).LastPoint;
				GAP.at(IDGap).active=false;
				GAP.at(IDGap).endTime=step;
				GapAux.IDList.insert(GapAux.IDList.end(), GAP.at(IDGap).IDList.begin(), GAP.at(IDGap).IDList.end());
				overlapGap=true;

				std::cout << "Gap overlap!! Middle point: "<< ID << std::endl;

			}

			ID=NormalizedID(GAP.at(IDGap).LastPoint+1, startIDRef);
			if ( initRef== ID)  // Se unen 
			{
				GapAux.InitPoint=GAP.at(IDGap).InitPoint;
				GAP.at(IDGap).active=false;
				GAP.at(IDGap).endTime=step;
				GapAux.IDList.insert(GapAux.IDList.end(), GAP.at(IDGap).IDList.begin(), GAP.at(IDGap).IDList.end());
				std::cout << "Gap overlap!! Middle point: "<< ID << std::endl;
				overlapGap=true;

			}
//		}
	  }	
		
		
	}
	return overlapGap;

}


void Monolayer::GetBoundListOK(std::vector<int> &ListAux , int ID )
{
	int IDCadh, IDAux,IDCadh2 ;

	bool state;
	ListAux.clear();

	if (!Nodes.at(ID)-> get_bound())
	{
		std::cout << "ID Causing trouble  "<<  ID <<std::endl;
		std::cout << "ListID Causing trouble:  " <<std::endl;
		for (int i=0; i <  ListAux.size(); i++)
			std::cout <<  ListAux[i] << ", ";
		std::cout <<std::endl;
	}

	IDCadh=Nodes.at(ID)-> CadhPointer.at(0);
//	IDCadh2=Nodes.at(ID)-> CadhPointer.at(1);
	state=Nodes.at(ID)-> get_bound();

	if (CadhPack.List.at(IDCadh).firstAdhPoint==ID )	// Check the point that the cadherin binds
		IDAux=CadhPack.List.at(IDCadh).secondAdhPoint;
	else
		IDAux=CadhPack.List.at(IDCadh).firstAdhPoint;

	ListAux.push_back(IDAux);

	if ( Nodes.at(ID)-> CadhPointer.size() > 1)
	{
		IDCadh=Nodes.at(ID)-> CadhPointer.at(1);   // Check if its bound to two points
		if (IDCadh > 0) // Tiene 2 unidos
		{
			if (CadhPack.List.at(IDCadh).firstAdhPoint==ID )
				IDAux=CadhPack.List.at(IDCadh).secondAdhPoint;
			else
				IDAux=CadhPack.List.at(IDCadh).firstAdhPoint;
	
			ListAux.push_back(IDAux);
	
		}
	}	
	
}

	 
/* //  Public-domain function by Darel Rex Finley, 2006.
double polygonArea(double *X, double *Y, int points) 
{

  double  area=0. ;
  int     i, j=points-1  ;

  for (i=0; i<points; i++) 
  {
    area+=(X[j]+X[i])*(Y[j]-Y[i]); 
    j=i; 
}

  return area*.5; 
}

*/

double Monolayer::polygonArea(std::vector<int>& pointsList) 
{

  double  area=0. ;
  int IDi, IDj;
  int     i;
  int    j=pointsList.size()-1  ;

 // std::set<int> myset;   // Si le metes duplicado te lo quita
 //std::vector<int> pointsListAux;

 // Remove duplicates
  std::vector<int>::iterator it;
  it = std::unique (pointsList.begin(), pointsList.end());   // 10 20 30 20 10 ?  ?  ?  ?
                                                         //                ^

  pointsList.resize( std::distance(pointsList.begin(),it) ); // 10 20 30 20 10

  for (i=0; i<pointsList.size(); i++) 
  {
  	IDi=pointsList[i];
  	IDj=pointsList[j];


    area+=(Nodes.at(IDj)->Coord[0]+Nodes.at(IDi)->Coord[0])*(Nodes.at(IDj)->Coord[1]- Nodes.at(IDi)->Coord[1]); 
    j=i; 
   }
   area=fabs(area*0.5);
  return area; 
}

double Monolayer::MinumunDistanceToVertex(std::vector<int> ListAux, int start, int end, int cellID, int& closestVertAux)
{
	int ID;
	double distAux=aristLength;;
	int distanceAdhPoint=aristLength;

	std::vector<int> VertAux;
	VertAux.push_back(GETVertexID(cellID, 0));
	VertAux.push_back(GETVertexID(cellID, 1));
	VertAux.push_back(GETVertexID(cellID, 2));
	VertAux.push_back(GETVertexID(cellID, 3));
	VertAux.push_back(GETVertexID(cellID, 4));
	VertAux.push_back(GETVertexID(cellID, 5));

	for (int i=0; i< ListAux.size(); i++)
	{
		ID=ListAux[i];
		for (int j=0; j< VertAux.size() ;j++)
		{
			if (VertAux[j]==ID)
				return double(0);
		}

	}

	for (int j=0; j< VertAux.size() ;j++)
	{
		ID=start;
		distAux=fabs(VertAux[j]-ID)*sepAdhPoint;
		if ( distAux < distanceAdhPoint )
		{
				distanceAdhPoint=distAux;
				closestVertAux=0;
		}

		ID=end;
		distAux=fabs(VertAux[j]-ID)*sepAdhPoint;
		if ( distAux < distanceAdhPoint )
		{
				distanceAdhPoint=distAux;
				closestVertAux=1;
		}

	}

	return distanceAdhPoint;
}

int Monolayer::NormalizedID(int ID, int start)
{
	if ( ID >= (start+nAdhPointCell))
		ID=ID-nAdhPointCell;
	else if (ID < start)
		ID=ID+nAdhPointCell;

	
	return ID;
}

void Monolayer::CalculateMembraneDistance()
{
	int cellID=floor(nCellRow/2);
	double celldistance, pointsDistance; 
	MembAdhData adhesionData;
	// Cell center
	int startID=GETVertexID(cellID,0);
	int endID=startID+nAdhPointCell;
	int ct=0;
	int auxCt;
	int auxID;

	for (auto PointID=startID; PointID < endID; ++PointID)
	{
		auxCt=0;
		for (int i=0; i < Pool.size(); ++i)
		{	
			if (i==cellID)  // Jump the cell in the middle. 
				i++;

			celldistance=NodesDist(PointID, Pool[i].GETIDCenter());
			if (celldistance < aristLength*2.2)
			{
				pointsDistance=MinimumMembraneDistance(PointID,i, auxID);
				adhesionData.distance=pointsDistance;
				adhesionData.refID=PointID;
				adhesionData.pointMinID=auxID;
				adhesionData.bound=Nodes.at(PointID) -> get_bound();
				auxCt++;					
			}
			if (auxCt==1) // Copy first case
				MembDistance[ct]=adhesionData;
			else if (MembDistance[ct].distance > adhesionData.distance) // In case distance is closer
				MembDistance[ct]=adhesionData;
		}
			
		ct++;

	}	
}

double Monolayer::MinimumMembraneDistance(int refPointID, int cellID, int &pointMin)  // returns Distance between refPoint and closest membrane point of cellID
{
	
	int auxCellPoint;
	double minimunDist, minimunDistAux;

	
	int startPoint=GETVertexID(cellID,0);
	pointMin=startPoint;
	minimunDist=NodesDist(refPointID, startPoint);
	for (auto i=1; i< nAdhPointCell; i++)
	{
		auxCellPoint=startPoint+i;
		minimunDistAux=NodesDist(refPointID, auxCellPoint);
		if (minimunDistAux < minimunDist)
		{
			minimunDist=minimunDistAux;
			pointMin=auxCellPoint;	
		}
	}
	return minimunDist;
}

double Monolayer::CalculateStress(Eigen::VectorXd &FInt, double EABalance, double lBalance , int a, int b)
{
	double forcemod, lCurrent, EACurrent;
	std::vector<double> VDir(2);

	lCurrent=NodesDist(a, b);
//	EACurrent=EABalance*lBalance/lCurrent;
//	forcemod=EACurrent*log(lCurrent/lBalance);
//	forcemod=EABalance/lBalance*( NodesDistX(a, b)-NodesDistInitX(a, b) );
	
//	forcemod=EABalance/lBalance*( lCurrent-lBalance);

	// Using EA as K
	forcemod=EABalance*( lCurrent-lBalance);


	return forcemod;
}



void Monolayer::RepForceCenter(Eigen::VectorXd &FInt)	// Force acting omn center by distance beetween centers
{
	double celldistance, forceRep;
	std::vector<double> VDir(2);
	int cellA, cellB;


	// Cell center
	for (int i=0; i < nCell-1; ++i)
	{
		for (int j=i+1; j < nCell; ++j)
		{
			cellA=Pool[i].GETIDCenter();
			cellB=Pool[j].GETIDCenter();
			celldistance=NodesDist(cellA, cellB);
			if (celldistance <  centerMinimumDist)
			{
				NodesVDir(cellA, cellB, VDir);
				forceRep=double(-1)*(celldistance-centerMinimumDist)*0.002;

				FInt(2*cellB)+=forceRep*VDir.at(0);
				FInt(2*cellB+1)+=forceRep*VDir.at(1);

				FInt(2*cellA)+=forceRep*VDir.at(0)*double(-1);
				FInt(2*cellA+1)+=forceRep*VDir.at(1)*double(-1);


			}
		}
	}

}

void Monolayer::RepulsionMembraneGlobalForce(Eigen::VectorXd &FInt)
{
	double distanceAux, KBendTotal, forceRep ;
	double auxRepulsion=double(1);
	int pA, pB;
	int cellIDAux, centerIDAux;
	std::vector<double> VDir(2);

	for(int i=0; i < MembDistanceGlobal.size(); i++)
	{
		pA=MembDistanceGlobal.at(i).firstPoint;
		pB=MembDistanceGlobal[i].secondPoint;
		distanceAux=NodesDist(pA, pB);
		MembDistanceGlobal[i].distance=distanceAux;
		cellIDAux=CellIDfromPoint(pA);
		centerIDAux=Pool.at(cellIDAux).GETIDCenter();

		if ((distanceAux >0.001) && (distanceAux < (CadhPack.GETlengthBal()/0.85)  ))
		{
			NodesVDir(pA, pB, VDir);
			KBendTotal=double(CadhPack.GETEA()*CadhPack.GETclusterMaxNumber()*0.55);
			forceRep=KBendTotal*(distanceAux-CadhPack.GETlengthBal())/mediumViscCoef*double(-1);

						// Add repulsion phenomenon
			
			auxRepulsion=CellOverlap(FInt,pA, pB, 1.);// In case one cell is inside another, change sign of the force
			forceRep=forceRep*auxRepulsion;

			FInt(2*pB)+=forceRep*VDir.at(0);
			FInt(2*pB+1)+=forceRep*VDir.at(1);

			FInt(2*pA)+=forceRep*VDir.at(0)*double(-1);
			FInt(2*pA+1)+=forceRep*VDir.at(1)*double(-1);
		}
		else  // If cell is inside another
		{
			if ( NodesDist(pA, centerIDAux) -2. > NodesDist(pB, centerIDAux)  )   // Center in closer to the membrane point of the other cell, apply repulsion
				CellOverlap(FInt,pA, pB, 2.1 );// In case one cell is inside another, change sign of the force

		}

	}
}


void Monolayer::RepulsionMembraneGlobalInteraction(Eigen::VectorXd &FInt)
{
	int IDaux, IDstart;
	std::vector<bool> MembDistGlobChecked;
	MembDistGlobChecked.resize((nAdhPointCell+1)*Pool.size(), false);

	MembDistanceGlobal.clear();

	for(int i=0; i < Pool.size(); i++)
	{
		IDstart=Pool[i].GETIDCenter()+1;
		for(int j=0; j < nAdhPointCell; j++)
		{
			IDaux=IDstart+j;
			if (! MembDistGlobChecked.at(IDaux)) // Check if it is not checked before
				CheckMembDistanceGlobal(IDaux, i, MembDistGlobChecked);
		}
	}
}

void Monolayer::CheckMembDistanceGlobal(int IDaux, int IDcell, std::vector<bool> &MembDistGlobChecked)
{
	double celldistance, distClosest, distClosestAux, distanceLimitCell;
	int IDclosest, IDclosestAux;
	
	 DistanceData MembDistanceGlobalAux;

	distanceLimitCell=aristLength*3;
	distClosest=distanceLimitCell*10;
	for (int i=0; i < nCell; ++i)
	{
	  if (IDcell!=i)
	  { 
		celldistance=NodesDist(Pool[i].GETIDCenter(), Pool[IDcell].GETIDCenter());
		if (celldistance < distanceLimitCell )
		{
			distClosestAux=CellClosestPoint(IDaux, i, IDclosestAux);
			if (distClosestAux < distClosest)
			{
				distClosest=distClosestAux;
				IDclosest=IDclosestAux;
			}

		}

	  }
	}

	MembDistanceGlobalAux.firstPoint=IDaux;
	MembDistanceGlobalAux.secondPoint=IDclosest;
	MembDistanceGlobalAux.distance=distClosest;
	MembDistanceGlobal.push_back(MembDistanceGlobalAux);

	MembDistGlobChecked.at(IDclosest)=true;  // We dont analyze it when it is reached

}

// Return distance to closest point in cell IDcellCheck, IDclosestAux 
double Monolayer::CellClosestPoint(int IDref, int IDcellCheck, int &IDclosestAux)
{
	int IDvertexAux, IDaux, IDstart;
	double distanceAux, distanceRef;
	distanceRef=10e10;

	 // Find the closest vertex
	for (int i=0; i < nArists; ++i)
	{
		distanceAux=NodesDist(IDref, GETVertexID(IDcellCheck,i));  // Distance between vertex
		if (distanceAux < distanceRef )
		{
			IDvertexAux=i;
			distanceRef=distanceAux;
		}
	}

	 // IDvertexAux
	distanceRef=10e10;
	IDstart=GETVertexID(IDcellCheck,IDvertexAux);


	// We have the closest vertex, now we go through both sides to see closest point 
	for (int i=1; i < nAdhPointSide; ++i)
	{
		IDaux=NormalizedID(IDstart+i, GETVertexID(IDcellCheck,0));
		distanceAux=NodesDist(IDref, IDaux);  // Distance between vertex
		if (distanceAux < distanceRef )
		{
			IDclosestAux=IDaux;
			distanceRef=distanceAux;
		}
	
	}

	for (int i=1; i < nAdhPointSide; ++i)
	{
		IDaux=NormalizedID(IDstart-i, GETVertexID(IDcellCheck,0));
		distanceAux=NodesDist(IDref, IDaux);  // Distance between vertex
		if (distanceAux < distanceRef )
		{
			IDclosestAux=IDaux;
			distanceRef=distanceAux;
		}
	
	}

	return 	distanceRef;
}





void Monolayer::RepulsionCellInCenter(Eigen::VectorXd &FInt)
{
	int auxIDCenter;
	int a;
	int b;
	int cellIDAux, centerIDAux;
	double lCurrent, forcemod, lBalance, KBendTotal;
	std::vector<double> VDir(2);
	
	for (int i=0; i < MembDistance.size(); i++)
	{
		a=MembDistance[i].refID;
		b=MembDistance[i].pointMinID;

		lCurrent=NodesDist(a, b);
	//	std::cout << "lCurrebt:  " << lCurrent << std::endl;

		if ((lCurrent >0.001) && (lCurrent < (CadhPack.GETlengthBal())  ))
		{


			lBalance=CadhPack.GETlengthBal();
			KBendTotal=double(CadhPack.GETEA()*CadhPack.GETclusterMaxNumber()*2.5);
			forcemod=KBendTotal*(lCurrent-lBalance)/mediumViscCoef;

			//repulsion applied in the radial direction
			cellIDAux=floor(a/(nAdhPointCell+1));
			centerIDAux=Pool.at(cellIDAux).GETIDCenter();
			NodesVDir(centerIDAux, a, VDir);
			
		//	std::cout << "Force:  " << forcemod << std::endl;
			FInt(2*a)+=forcemod*VDir.at(0);
			FInt(2*a+1)+=forcemod*VDir.at(1);

			//repulsion applied in the radial direction
			cellIDAux=floor(b/(nAdhPointCell+1));
			centerIDAux=Pool.at(cellIDAux).GETIDCenter();
			NodesVDir(centerIDAux,b, VDir);
			
		//	std::cout << "Force:  " << forcemod << std::endl;
			FInt(2*b)+=forcemod*VDir.at(0);
			FInt(2*b+1)+=forcemod*VDir.at(1);

		}
	}

}

void Monolayer::CalculateNodesForce(Eigen::VectorXd &FInt, double EABalance, double barVisc, double lBalance , int a, int b)
{   // a es boundry
	double forcemod, lCurrent, EACurrent;
	std::vector<double> VDir(2);

	lCurrent=NodesDist(a, b);
	EACurrent=EABalance*lBalance/lCurrent;
 //	forcemod=EACurrent*log(lCurrent/lBalance);
 //	forcemod=EABalance/lBalance*( NodesDistX(a, b)-NodesDistInitX(a, b) );
	
 //	forcemod=EABalance/lBalance*( lCurrent-lBalance);

	// Using EA as K
	forcemod=EABalance*(lCurrent-lBalance)/barVisc;


	NodesVDir(a,b, VDir);

 /*	// For small deformation
	double c=NodesDistInitX(a, b)/NodesDistInit(a, b);
	double s=NodesDistInitY(a, b)/NodesDistInit(a, b);
	NodesInitVDir(a,b, VDir);   // b-a
 */
 //	forcemod=EA/lBalance*( NodesDistX(a, b)-NodesDistInitX(a, b) );
	FInt(2*b)+=forcemod*VDir.at(0)*double(-1);
 //	forcemod=EA/lBalance*( NodesDistY(a, b)-NodesDistInitY(a ,b) );
	FInt(2*b+1)+=forcemod*VDir.at(1)*double(-1);

	//	forcemod=EA/lBalance*( NodesDistX(a, b)-NodesDistInitX(a, b) );
		FInt(2*a)+=forcemod*VDir.at(0);
	//	forcemod=EA/lBalance*( NodesDistY(a, b)-NodesDistInitY(a, b) );
		FInt(2*a+1)+=forcemod*VDir.at(1);
	
}




void Monolayer::ContractionByRepulsion(Eigen::VectorXd &FInt, int ID, double multiplier)
{
	int centerIDAux, cellIDAux, startIDRef;

	cellIDAux=CellIDfromPoint(ID);
	centerIDAux=Pool.at(cellIDAux).GETIDCenter();

	std::vector<double> VDir(2);
	NodesVDir(ID, centerIDAux, VDir);
	FInt(2*ID)+=VDir.at(0)*maxForceRep*0.3*multiplier/(mediumViscCoef+bdryVisc);
	FInt(2*ID+1)+=VDir.at(1)*maxForceRep*0.3*multiplier/(mediumViscCoef+bdryVisc) ;
}

double Monolayer::CellOverlap(Eigen::VectorXd &FInt, int a, int b, double multiplier)
{
	int centerIDAux, cellIDAux, ID, startIDRef;

	cellIDAux=CellIDfromPoint(a);
	startIDRef=GETVertexID(cellIDAux,0);
	std::vector<double> VDir(2);
	centerIDAux=Pool.at(cellIDAux).GETIDCenter();
	//  Remove comment for contraction as repulsion
		ContractionByRepulsion(FInt,a, multiplier);
		for (int i=1; i <= 1;i++)
		{
			ID=NormalizedID(a-i, startIDRef);
			ContractionByRepulsion(FInt,ID, multiplier);
			ID=NormalizedID(a+i, startIDRef);
			ContractionByRepulsion(FInt,ID, multiplier);
		}
		
		cellIDAux=CellIDfromPoint(b);
		startIDRef=GETVertexID(cellIDAux,0);

		ContractionByRepulsion(FInt,b, multiplier);
		for (int i=1; i <= 1;i++)
		{
			ID=NormalizedID(b-i, startIDRef);
			ContractionByRepulsion(FInt,ID, multiplier);
			ID=NormalizedID(b+i, startIDRef);
			ContractionByRepulsion(FInt,ID, multiplier);
		}
	

	if ( NodesDist(a, centerIDAux) -10. > NodesDist(b, centerIDAux)  )   // Center in closer to the membrane point of the other cell, apply repulsion
	{
		// Contract both arms
/*		//  Remove comment for contraction as repulsion
		ContractionByRepulsion(FInt,a);
		for (int i=1; i <= 1;i++)
		{
			ID=NormalizedID(a-i, startIDRef);
			ContractionByRepulsion(FInt,ID);
			ID=NormalizedID(a+i, startIDRef);
			ContractionByRepulsion(FInt,ID);
		}
		
		cellIDAux=CellIDfromPoint(b);
		startIDRef=GETVertexID(cellIDAux,0);

		ContractionByRepulsion(FInt,b);
		for (int i=1; i <= 1;i++)
		{
			ID=NormalizedID(b-i, startIDRef);
			ContractionByRepulsion(FInt,ID);
			ID=NormalizedID(b+i, startIDRef);
			ContractionByRepulsion(FInt,ID);
		}
	*/
		return double(-1.);
	}
	else 
		return double(1);

}


void Monolayer::CalculateNodesForceCadh(Eigen::VectorXd &FInt, double EABalance, double barVisc, double lBalance , int a, int b)
{   // a es boundry
	double forcemod, lCurrent, EACurrent;
	double auxRepulsion=double(1);  // In case one cell is inside another, change sign of the force
	std::vector<double> VDir(2);

	lCurrent=NodesDist(a, b);
	EACurrent=EABalance*lBalance/lCurrent;
 
	// Using EA as K
	forcemod=EABalance*(lCurrent-lBalance)/barVisc;

/*
	// Add repulsion phenomenon
	if (forcemod < 0 )
	{
		auxRepulsion=CellOverlap(FInt,a, b);// In case one cell is inside another, change sign of the force
		forcemod=forcemod*3.*auxRepulsion;
	}
*/
	NodesVDir(a,b, VDir);

	FInt(2*b)+=forcemod*VDir.at(0)*double(-1);
 
	FInt(2*b+1)+=forcemod*VDir.at(1)*double(-1);
	
	FInt(2*a)+=forcemod*VDir.at(0);

	FInt(2*a+1)+=forcemod*VDir.at(1);
	
}
void Monolayer::CalculateNodesForceBdry(Eigen::VectorXd &FInt, double EABalance, double barVisc, double lBalance , int a, int b)
{   // a es boundry
	double forcemod, lCurrent, EACurrent;
	std::vector<double> VDir(2);

	lCurrent=NodesDist(a, b);
	EACurrent=EABalance*lBalance/lCurrent;
 //	forcemod=EACurrent*log(lCurrent/lBalance);
 //	forcemod=EABalance/lBalance*( NodesDistX(a, b)-NodesDistInitX(a, b) );
	
 //	forcemod=EABalance/lBalance*( lCurrent-lBalance);

	// Using EA as K
	forcemod=EABalance*(lCurrent-lBalance);

	
	NodesVDir(a,b, VDir);

 /*	// For small deformation
	double c=NodesDistInitX(a, b)/NodesDistInit(a, b);
	double s=NodesDistInitY(a, b)/NodesDistInit(a, b);
	NodesInitVDir(a,b, VDir);   // b-a
 */
 //	forcemod=EA/lBalance*( NodesDistX(a, b)-NodesDistInitX(a, b) );
	FInt(2*b)+=forcemod*VDir.at(0)*double(-1)/barVisc;
 //	forcemod=EA/lBalance*( NodesDistY(a, b)-NodesDistInitY(a ,b) );
	FInt(2*b+1)+=forcemod*VDir.at(1)*double(-1)/barVisc;

//	if(!bdry)
//	{
	//	forcemod=EA/lBalance*( NodesDistX(a, b)-NodesDistInitX(a, b) );
		FInt(2*a)+=forcemod*VDir.at(0)/(barVisc + nucleusViscCoef );   // Adding more viscosity to nucleus movememtn
	//	forcemod=EA/lBalance*( NodesDistY(a, b)-NodesDistInitY(a, b) );
		FInt(2*a+1)+=forcemod*VDir.at(1)/(barVisc + nucleusViscCoef);
//	}



}


void Monolayer::CalculateForceInternAlternative(SpMat &KMat, Eigen::VectorXd &FInt, Eigen::VectorXd &XIncTotal )
{
	FInt=KMat*XIncTotal*double(-1);
}



void Monolayer::UpdateNodesPosition(Eigen::VectorXd &XIncEfeective)
{	
	for (int i=0; i<meshSize; ++i)
	{	
		Nodes.at(i)->Coord.at(0)=Nodes.at(i)->Coord.at(0)+XIncEfeective(2*i);
		Nodes.at(i)->Coord.at(1)=Nodes.at(i)->Coord.at(1)+XIncEfeective(2*i+1);
		//(it->second)->Coord.at(0)=NodesInitialPos.at(i).at(0)+XIncEfeective(2*i);
		//(it->second)->Coord.at(1)=NodesInitialPos.at(i).at(0)+XIncEfeective(2*i+1);
	}
}

void Monolayer::UpdateNodesPosition(std::vector<double> &XIncEfeective)
{	
	for (int i=0; i<meshSize; ++i)
	{
		Nodes.at(i)->Coord.at(0)=Nodes.at(i)->Coord.at(0)+XIncEfeective[2*i];
		Nodes.at(i)->Coord.at(1)=Nodes.at(i)->Coord.at(1)+XIncEfeective[2*i+1];
		//(it->second)->Coord.at(0)=NodesInitialPos.at(i).at(0)+XIncEfeective(2*i);
		//(it->second)->Coord.at(1)=NodesInitialPos.at(i).at(0)+XIncEfeective(2*i+1);
	}
}

void Monolayer::CheckNodesIncrement()
{
	
	for (int i=0; i<meshSize; ++i)
	{
		
		XIncTotalCheck(2*i)=Nodes.at(i)->Coord.at(0)-NodesInitialPos.at(i).at(0);
		XIncTotalCheck(2*i+1)=Nodes.at(i)->Coord.at(1)-NodesInitialPos.at(i).at(1);
	}

}

/*
SpMat BuildGlobalK()   // Change with new way of storing Sparse matrix
{

	SpMat MatrixKFixed(2*meshSize,2*meshSize);

	std::vector<T> CadherinTriplet, tripletList;

	CadherinTriplet=BuildTripletCadhConnectivity();

	tripletList.reserver(CadherinTriplet.size()+ tripletListMembBdry.size());

	for (int i=0; i< tripletListMembBdry.size();++i)
		tripletList.push_back(tripletListMembBdry.at(i));

	for (int i=0; i< CadherinTriplet.size();++i)
		tripletList.push_back(CadherinTriplet.at(i));

	MatrixKFixed.setFromTriplets(tripletList.begin(), tripletList.end());
	return MatrixKFixed;
}
*/

void Monolayer::ForceGeneration()
{
	bool radialF=true;
	double rndRad, rndCt, rndFlow, rndProt, auxDir, rndBr;
	Eigen::VectorXd ForceStepVc(2*meshSize);
	std::vector<double> vDir(2);
	std::uniform_real_distribution<double> distribution(0, 1);
	int ctPropMemb=0,ctPropBdry=0, ctPropProt=0;

	if (forceCt >= stepForce)
	{
		forceCt=1;

		// For radial
		ForceStepVc= Eigen::VectorXd::Zero(2*meshSize);
	
		for (int i=0; i<ConnectivityBdry.at(0).size()-extraBdry;++i) // Remember (0 -> bdry; 1-> Membpoint)
		{

			if (ctPropBdry==0 || ctPropBdry==nPropagationBdry )
			{ 
				ctPropBdry=0;
				rndRad=distribution(generator)*maxRadForce*ratioZeroRad;
				if (rndRad > maxRadForce )
					rndRad=distribution(generator)*maxRadForce*0.1;

				// remove for normal working. With this we can only have maximun or 0 force			
				else
					rndRad=maxRadForce;

				rndRad=rndRad/(mediumViscCoef+bdryVisc);

			}	
			ctPropBdry++;

			NodesVDir(ConnectivityBdry.at(1).at(i), ConnectivityBdry.at(0).at(i), vDir);
			if( nArists==1  && (!cellCenter) )   // In case of wall config not consider bdry for y axis when apply contractile force
			{
				ForceStepVc(ConnectivityBdry.at(1).at(i)*2)=vDir.at(0)*rndRad; // x direction
			}
			else
			{
				ForceStepVc(ConnectivityBdry.at(1).at(i)*2)=vDir.at(0)*rndRad; // x direction
				ForceStepVc(ConnectivityBdry.at(1).at(i)*2+1)=vDir.at(1)*rndRad; // y direction	
			}

			if (testFormulation)
			{
				if (step==9)
					rndRad=double(0);
				ForceStepVc(ConnectivityBdry.at(1).at(i)*2)=rndRad; // x direction
				rndRad=distribution(generator);
				if (step==9)
					rndRad=double(0);
				ForceStepVc(ConnectivityBdry.at(1).at(i)*2+1)=0; // y direction	
			}
	
		//	ForceInc(ConnectivityBdry.at(1).at(i)*2)=ForceStepVc(ConnectivityBdry.at(1).at(i)*2)
		//											 -ForceLastStep (ConnectivityBdry.at(1).at(i)*2);
			
		}
		// For Protrusion
	
		for (int i=0; i<ConnectivityBdry.at(0).size()-extraBdry;++i) // Remember (0 -> bdry; 1-> Membpoint)
		{

			if (ctPropProt==0 || ctPropProt==nPropagationProt )
			{ 
				ctPropProt=0;
				rndProt=distribution(generator)*maxProtForce*ratioZeroProt;
				if (rndProt > maxProtForce)
					rndProt=double(0);

				// remove for normal working. With this we can only have maximun or 0 force			
				else
					rndProt=maxProtForce;

				rndProt=rndProt/(mediumViscCoef+bdryVisc)*double(-1);

			}	
			ctPropProt++;

			NodesVDir(ConnectivityBdry.at(1).at(i), ConnectivityBdry.at(0).at(i), vDir);
			if( nArists==1  && (!cellCenter) )   // In case of wall config not consider bdry for y axis when apply contractile force
			{
				ForceStepVc(ConnectivityBdry.at(1).at(i)*2)+=vDir.at(0)*rndProt; // x direction
			}
			else
			{
				ForceStepVc(ConnectivityBdry.at(1).at(i)*2)+=vDir.at(0)*rndProt; // x direction
				ForceStepVc(ConnectivityBdry.at(1).at(i)*2+1)+=vDir.at(1)*rndProt; // y direction	
			}
			
		}

		

		// For cortical
		for (int i=0; i<ConnectivityMemb.at(0).size();++i) // 
		{
			

			if (ctPropMemb==0 || ctPropMemb==nPropagationMemb )
			{
				ctPropMemb=0;
				rndCt=distribution(generator)*maxCtForce*ratioZeroCt;
				if (rndCt > maxCtForce )
					rndCt=double(0);
				// remove for normal working. With this we can only have maximun or 0 force			
				else
					rndRad=maxCtForce;


				rndCt=rndCt/(mediumViscCoef+membVisc );

			}	
			ctPropMemb++;


			NodesVDir(ConnectivityMemb.at(0).at(i), ConnectivityMemb.at(1).at(i), vDir);
			if( i==0 )   // In case of wall config not consider bdry for y axis when apply contractile force
			{
				ForceStepVc(ConnectivityMemb.at(0).at(i)*2)+=vDir.at(0)*rndCt; // x direction
				ForceStepVc(ConnectivityMemb.at(0).at(i)*2+1)+=vDir.at(1)*rndCt; // y direction	
			}
			
			ForceStepVc(ConnectivityMemb.at(1).at(i)*2)+=vDir.at(0)*rndCt; // x direction
			ForceStepVc(ConnectivityMemb.at(1).at(i)*2+1)+=vDir.at(1)*rndCt; // y direction	
		}
	
		// For directional
		dirFlowForce=dirFlowForce*3.1416/180;
		vDir.at(0)=cos(dirFlowForce);
		vDir.at(1)=sin(dirFlowForce);
		for (int i=0; i<ConnectivityBdry.at(0).size();++i) // Remember (0 -> bdry; 1-> Membpoint)
		{
			rndFlow=distribution(generator)*maxFlowForce;
		/*	
			if( i==0 )   // In case of wall config not consider bdry for y axis when apply contractile force
			{
				ForceStepVc(ConnectivityBdry.at(0).at(i)*2)+=vDir.at(0)*rndF; // x direction
				ForceStepVc(ConnectivityBdry.at(0).at(i)*2+1)+=vDir.at(1)*rndF; // y direction	
			}
		*/	
			rndFlow=rndFlow/mediumViscCoef ;
			ForceStepVc(ConnectivityBdry.at(1).at(i)*2)+=vDir.at(0)*rndFlow; // x direction
			ForceStepVc(ConnectivityBdry.at(1).at(i)*2+1)+=vDir.at(1)*rndFlow; // y direction	
		}
	
		ForceInc=(ForceStepVc-ForceLastStep)/stepSmoothForce;
		ForceStep=ForceStep+ForceInc;
	//	print_data_intern2(ForceStepVc);
		ForceLastStep=ForceStepVc;
	}

	else
	{
		if (forceCt < stepSmoothForce )
			ForceStep=ForceStep+ForceInc;
		// Else ForceStep does not change any more


		forceCt++;
	}

/*
	// For brownian dynamics
	if (stepBrownian == 10 )
	{	
		for (int i=0; i<ConnectivityBdry.at(0).size()-extraBdry;++i) // Remember (0 -> bdry; 1-> Membpoint)
		{
	
		//	rndBr=distribution(generator)*(maxRadForce+maxCtForce)*0.05;
			rndBr=distribution(generator)*0.105;
			rndBr=rndBr/(mediumViscCoef+bdryVisc);
			auxDir=distribution(generator)*2*PI;
			vDir.at(0)=cos(auxDir);
			vDir.at(1)=sin(auxDir);		
		
			ForceStepBr(ConnectivityBdry.at(1).at(i)*2)=vDir.at(0)*rndBr; // x direction
			ForceStepBr(ConnectivityBdry.at(1).at(i)*2+1)=vDir.at(1)*rndBr; // y direction	
		
		}
		stepBrownian=0;
	}
	// 
	else
		stepBrownian++;

*/


	//   MakeUnboundForceZero();   // Uncomment to make force=0 in unbound nodes 
}

// Function to make force=0 in unbound nodes 
// Function to make force=0 in unbound nodes 
void Monolayer::MakeUnboundForceZero()   // Optimize with unbound
{
	
	int nodeIDAux; 
	for (int i=0; i<ConnectivityMemb.at(1).size();++i)  // Remember (0 -> bdry; 1-> Membpoint)
	{
		nodeIDAux=ConnectivityMemb[1][i];
		if (Nodes[nodeIDAux]-> get_bound()==false)
		{
			ForceStep(nodeIDAux*2)=0.; // x direction
			ForceStep(nodeIDAux*2+1)=0.; // y direction	
		}
			
		
	}
}




void Monolayer::BuildSparseMatrixMembBdry() // Used in constructor
{
	int i,a,b;
	double effectiveCoef;
	std::vector<T> KTriplet;
	std::vector<T> CTriplet;
	KMemBdry.resize(2*meshSize,2*meshSize);
	CMemBdry.resize(2*meshSize,2*meshSize);

	for (i=0; i< ConnectivityMemb.at(0).size();++i)
	{
		a=ConnectivityMemb.at(0).at(i);
		b=ConnectivityMemb.at(1).at(i);
		TripletXY(KTriplet, CTriplet, membEA, membVisc, NodesDistInit(a,b),	a, b, false);
	}
	

	for (i=0; i< ConnectivityBdry.at(0).size();++i)
	{
		a=ConnectivityBdry.at(0).at(i);
		b=ConnectivityBdry.at(1).at(i);
		TripletXY(KTriplet, CTriplet, bdryEA, bdryVisc, NodesDistInit(a,b),	a, b, true);
	}

	KMemBdry.setFromTriplets(KTriplet.begin(), KTriplet.end());
	CMemBdry.setFromTriplets(CTriplet.begin(), CTriplet.end());
}


void Monolayer::BuilSparseMatrixCadhConnectivity(SpMat &KMat, SpMat &CMat)
{
	double EA;
	std::vector<T> KTriplet;
	std::vector<T> CTriplet;
	KMat.resize(2*meshSize,2*meshSize);
	CMat.resize(2*meshSize,2*meshSize);
	int i, CadhID;
	for (i=0; i<CadhPack.ActiveCadh.size();++i)
	{
		CadhID=CadhPack.ActiveCadh.at(i);
	
		EA=CadhPack.GETEA()*CadhPack.List.at(CadhID).cadherinNumber;
		TripletXY(KTriplet, CTriplet, EA, CadhPack.GETvisc(), 
				  CadhPack.GETlengthBal(), CadhPack.List.at(CadhID).firstAdhPoint, 
						CadhPack.List.at(CadhID).secondAdhPoint, false );

	}
	
	KMat.setFromTriplets(KTriplet.begin(), KTriplet.end());
	CMat.setFromTriplets(CTriplet.begin(), CTriplet.end());
}

void Monolayer::TripletXY(std::vector<T>& KTriplet,std::vector<T>& CTriplet,  double EA, 
							double visc, double restLength, int a, int b, bool bdry)
{
	 double effectiveCoef, s,c ;

	 if (smallDefResolution==1)
	 {
	 	c=NodesDistInitX(a, b)/NodesDistInit(a, b);
	 	s=NodesDistInitY(a, b)/NodesDistInit(a, b);
	 }
	 else //Update matrix depending on current position
	 {
	 	c=NodesDistX(a, b)/NodesDist(a, b);
	 	s=NodesDistY(a, b)/NodesDist(a, b);
	 }

	 if (bdry)
	 {
	 	//effectiveCoef=EA/restLength;
	 	effectiveCoef=KMatrixCoef(NodesDist(a, b),restLength,EA );

	 	BuildTripletLocalKBdry(KTriplet,  2*b, effectiveCoef, c,s );

	 	effectiveCoef=visc;
	 	BuildTripletLocalKBdry(CTriplet,  2*b, effectiveCoef, c,s );
	 }

	 else

	 {
	 	//effectiveCoef=EA/restLength;
	 	effectiveCoef=KMatrixCoef(NodesDist(a, b),restLength,EA );
	 	BuildTripletLocalK(KTriplet, 2*a , 2*b, effectiveCoef, c,s );

	 	effectiveCoef=visc;
	 	BuildTripletLocalK(CTriplet, 2*a , 2*b, effectiveCoef, c,s );
	 }

}

double Monolayer::KMatrixCoef(double lCurrent,double lBalance, double EABalance)
{
	double Kcoef, EACurrent, logCoef;

//	 if (smallDefResolution==1)
	 if (smallDefResolution<2)
	 {
	 	Kcoef=EABalance/lBalance;
	 }
	 else
	 {
	/*
			sigma=E.at(i)*log(L_current/L_initial);
			a=L_initial*A/L_current;    // current area of the bar
			K_l=a*sigma/L_current + a/L_current*(E.at(i)-2*sigma);
	*/
		logCoef=log(lCurrent/lBalance);
		EACurrent=EABalance*lBalance/lCurrent;
		Kcoef=EACurrent*logCoef/lCurrent + EACurrent/lCurrent*(1-2*logCoef);
	}
	return Kcoef;
}

/*
void Monolayer::TripletXYCadh(std::vector<T>& KTriplet,std::vector<T>& CTriplet,  int a, int b)
{
	double effectiveCoef;
	// X positiom
	 effectiveCoef=CadhPack.GETEA()/CadhPack.GETlengthBal()*NodesDistX(a, b)/NodesDist(a, b);
	 BuildTripletLocalK(KTriplet, 2*a , 
							2*b, effectiveCoef );
	effectiveCoef=CadhPack.GETvisc()*NodesDistX(a, b)/NodesDist(a, b); 
	 BuildTripletLocalK(CTriplet, 2*a , 
							2*b, effectiveCoef ); 
	// Y position, visc
	effectiveCoef=CadhPack.GETEA()/CadhPack.GETlengthBal()*NodesDistY(a, b)/NodesDist(a, b);
	  BuildTripletLocalK(KTriplet, 2*a+1 , 
							2*b+1, effectiveCoef );
	effectiveCoef=CadhPack.GETvisc()*NodesDistY(a, b)/NodesDist(a, b); 
	 BuildTripletLocalK(CTriplet, 2*a+1 , 
							2*b+1, effectiveCoef ); 
}

void Monolayer::TripletXYMemb(std::vector<T>& KTriplet,std::vector<T>& CTriplet, int a, int b)
{
	double effectiveCoef;
	// X positiom
	 effectiveCoef=membEA/NodesDistInit(a,b)*NodesDistX(a, b)/NodesDist(a, b);
	 BuildTripletLocalK(KTriplet, 2*a , 
							2*b, effectiveCoef );
	effectiveCoef=membVisc*NodesDistX(a, b)/NodesDist(a, b); 
	BuildTripletLocalK(CTriplet, 2*a , 
							2*b, effectiveCoef ); 
	// Y position, visc
	effectiveCoef=membEA/NodesDistInit(a,b)*NodesDistY(a, b)/NodesDist(a, b);
	 BuildTripletLocalK(KTriplet, 2*a+1 , 
							2*b+1, effectiveCoef );
	effectiveCoef=membVisc*NodesDistY(a, b)/NodesDist(a, b); 
	BuildTripletLocalK(CTriplet, 2*a+1 , 
							2*b+1, effectiveCoef ); 
}

void Monolayer::TripletXYBdry(std::vector<T>& KTriplet,std::vector<T>& CTriplet,  int a, int b)
{
	double effectiveCoef;
	// X positiom
	 effectiveCoef=bdryEA/NodesDistInit(a,b)*NodesDistX(a, b)/NodesDist(a, b);
	 KTriplet.push_back(T( 2*b , 
							2*b, effectiveCoef ));
	effectiveCoef=bdryVisc*NodesDistX(a, b)/NodesDist(a, b); 
	CTriplet.push_back(T( 2*b , 
							2*b, effectiveCoef )); 
	// Y position, visc
	effectiveCoef=bdryEA/NodesDistInit(a,b)*NodesDistY(a, b)/NodesDist(a, b);
	 KTriplet.push_back(T( 2*b+1 , 
							2*b+1, effectiveCoef ));
	effectiveCoef=bdryVisc*NodesDistY(a, b)/NodesDist(a, b); 
	CTriplet.push_back(T( 2*b+1 , 
							2*b+1, effectiveCoef )); 
}

*/

void Monolayer::BuildTripletLocalK( std::vector<T>& tripletList, int a, int b, double coeff, double c, double s)
{

	tripletList.push_back(T(a , a, coeff*c*c ));  
	tripletList.push_back(T(a+1 , a+1, coeff*s*s ));  
	tripletList.push_back(T(a+1 , a, coeff*c*s ));  
	tripletList.push_back(T(a , a+1, coeff*c*s ));

	tripletList.push_back(T(b , b, coeff*c*c ));  
	tripletList.push_back(T(b+1 , b+1, coeff*s*s ));  
	tripletList.push_back(T(b+1 , b, coeff*c*s ));  
	tripletList.push_back(T(b , b+1, coeff*c*s )); 

	tripletList.push_back(T(a , b, -coeff*c*c ));  
	tripletList.push_back(T(a+1 , b+1, -coeff*s*s ));  
	tripletList.push_back(T(a+1 , b, -coeff*c*s ));  
	tripletList.push_back(T(a , b+1, -coeff*c*s ));  

	tripletList.push_back(T(b , a, -coeff*c*c ));  
	tripletList.push_back(T(b+1 , a+1, -coeff*s*s ));  
	tripletList.push_back(T(b+1 , a, -coeff*c*s ));  
	tripletList.push_back(T(b , a+1, -coeff*c*s )); 

}
void Monolayer::BuildTripletLocalKBdry( std::vector<T>& tripletList, int b, double coeff, double c, double s)
{

	tripletList.push_back(T(b , b, coeff*c*c ));  
	tripletList.push_back(T(b+1 , b+1, coeff*s*s ));  
	tripletList.push_back(T(b+1 , b, coeff*c*s ));  
	tripletList.push_back(T(b , b+1, coeff*c*s )); 

}

void Monolayer::VDirCalculation(std::vector<double> a, std::vector<double> b, 
								std::vector<double> &VDir)
{

}

double Monolayer::Vmod(std::vector<double> a, std::vector<double> b)
{


}

bool Monolayer::IsMajorInCell(int a, int b) // Return true if a is major or equal
{
	if ( (a-b) == 0  )
		return true;
	else if( (a-b) > 0  ) 
	{
		if ( fabs (a-b) < floor(nAdhPointCell/6) )
			return true;
		else 
			return false;
	}
	else // (a-b) < 0
	{
		if ( fabs (a-b) < floor(nAdhPointCell/6) )
			return false;
		else 
			return true;
	}
} 

int Monolayer::ReturnMajorInCell(int a, int b) // Return major point
{
	if( (a-b) >= 0  )
	{
		if ( fabs (a-b) < floor(nAdhPointCell/6) )
			return a;
		else 
			return b;
	}
	else // (a-b) < 0
	{
		if ( fabs (a-b) < floor(nAdhPointCell/6) )
			return b;
		else 
			return a;
	}	
} 




double Monolayer::NodesDist(int aIndex, int bIndex)  // b-a
{
	double length=double(0);
	double xa, xb;	
	for (int i=0; i<(Nodes.at(aIndex))-> Coord.size();++i)
	{
		xb=(Nodes.at(bIndex))-> Coord.at(i);
		xa=(Nodes.at(aIndex))-> Coord.at(i);

		length+=pow( xb-xa, 2 ); 
	}

	length=sqrt(length);
	return length;
}
double Monolayer::NodesDistInit(int aIndex, int bIndex)  // b-a
{
	double length=double(0);
	double xa, xb;
	for (int i=0; i<NodesInitialPos.at(bIndex).size();++i)
	{
		xb=NodesInitialPos.at(bIndex).at(i);
		xa=NodesInitialPos.at(aIndex).at(i);

		length+=pow( xb-xa, 2 ); 
	}
	length=sqrt(length);
	return length;
}

double Monolayer::NodesAngleInit(int aIndex, int bIndex)  // b-a
{
	double angleVal=double(0);
	double distX, distY;
	
	distX=NodesDistInitX(aIndex, bIndex); // b- a
	distY=NodesDistInitY(aIndex, bIndex);
	angleVal=atan2(distY,distX);
	
	return angleVal;
}



double Monolayer::NodesAngle(int aIndex, int bIndex)  // b-a
{
	double angleVal=double(0);
	double distX, distY;
	
	distX=NodesDistX(aIndex, bIndex); // b- a
	distY=NodesDistY(aIndex, bIndex);
	angleVal=atan2(distY,distX);
	
	return angleVal;
}

double Monolayer::NodesDistX(int aIndex, int bIndex)  // b-a
{
	double length;
		length= (Nodes.at(bIndex))-> Coord.at(0)-(Nodes.at(aIndex))-> Coord.at(0); 
//	length=fabs(length);
	return length;
}
double Monolayer::NodesDistY(int aIndex, int bIndex)  // b-a
{
	double length;
		length= (Nodes.at(bIndex))-> Coord.at(1)-(Nodes.at(aIndex))-> Coord.at(1); 
//	length=fabs(length);
	return length;
}
double Monolayer::NodesDistInitX(int aIndex, int bIndex)  // b-a
{
	double length;
		length= NodesInitialPos.at(bIndex).at(0)-NodesInitialPos.at(aIndex).at(0); 
//	length=fabs(length);
	return length;
}
double Monolayer::NodesDistInitY(int aIndex, int bIndex)  // b-a
{
	double length;
		length= NodesInitialPos.at(bIndex).at(1)-NodesInitialPos.at(aIndex).at(1); 
//	length=fabs(length);
	return length;
}

void Monolayer::AngleVDir(double angle, std::vector<double> &VDir) // b-a
{
	VDir[0]=cos(angle);
	VDir[1]=sin(angle);
}

void Monolayer::NodesVDir(int aIndex, int bIndex, std::vector<double> &VDir) // b-a
{
	double length=double(0);
	for (int i=0; i<(Nodes.at(aIndex))-> Coord.size();++i)
		length+=pow( (Nodes.at(bIndex))-> Coord.at(i)-(Nodes.at(aIndex))-> Coord.at(i), 2); 
	length=sqrt(length);
	for (int i=0; i<(Nodes.at(aIndex))-> Coord.size();++i)
		VDir.at(i)=( (Nodes.at(bIndex))-> Coord.at(i)-(Nodes.at(aIndex))-> Coord.at(i) )/length;
}
void Monolayer::NodesInitVDir(int aIndex, int bIndex, std::vector<double> &VDir) // b-a
{
	double length=double(0);
	for (int i=0; i<NodesInitialPos.at(aIndex).size();++i)
		length+=pow( NodesInitialPos.at(bIndex).at(i)-NodesInitialPos.at(aIndex).at(i), 2); 
	length=sqrt(length);
	for (int i=0; i<NodesInitialPos.at(aIndex).size();++i)
		VDir.at(i)=( NodesInitialPos.at(bIndex).at(i)-NodesInitialPos.at(aIndex).at(i) )/length;
}

void Monolayer::CalculateMiddlePoint(int aIndex, int bIndex, std::vector<double> &middlePoint) // b-a
{
	std::vector<double> VDir;
	for (int i=0; i<(Nodes.at(aIndex))-> Coord.size();++i)
		VDir.at(i)= NodesInitialPos.at(bIndex).at(i)-NodesInitialPos.at(aIndex).at(i); 
	
	for (int i=0; i<(Nodes.at(aIndex))-> Coord.size();++i)
		middlePoint.at(i)= (Nodes.at(aIndex))-> Coord.at(i) + VDir.at(i)/double(2);
}


int Monolayer::GETVertexID(int cellID, int vertex)
{
	int returnID;
	returnID=cellID*(nAdhPointCell+nCenterPointCell)+1;   // Placed in vertex number 0
	returnID=returnID+vertex*nAdhPointSide;
	return returnID;
}


/*
void Monolayer::UpdateCadherins()
{
	 for(auto it=CadhPack.List.begin(); it!=CadhPack.List.end(); ++it )
	 {
	 	(*it).length= NodesDist( (*it).firstAdhPoint, (*it).SecondAdhPoint );
	 	//(*it).force=(*it).length*CadhPack.GETK()/CadhPack.GETlengthBal();
	 }
	
}
*/
void Monolayer::printDataStep(char *index_cluster)
{

	char name[100];
	char orden[1024]; 
	char * extension=".txt";
	FILE*fOut;
	
	sprintf(name, "MembraneDistance_%s%s" , index_cluster, extension);				// 
	fOut = fopen(name,"a");   //para reescribir se usa "a"

	for (auto it=CadhPack.List.begin(); it !=CadhPack.List.end(); )
	{
		fprintf (fOut, "%4.1f ", (*it).length);  		    
		it++;
	}

	fprintf(fOut, "\n");
	fclose(fOut);

	sprintf(name, "GapTotalTime_%s%s" , index_cluster, extension);
	fOut = fopen(name,"a");   //para reescribir se usa "a"

	for (auto it=CadhPack.List.begin(); it !=CadhPack.List.end(); )
	{
		fprintf (fOut, "%4.1f ", (*it).gapTotTime);   // Print mean traction and last traction		    
		it++;
	}

	fprintf(fOut, "\n");
	fclose(fOut);

}

void Monolayer::printDataStepVertex(char *index_cluster,int t)
{

	char name[100];
	char orden[1024]; 
	char * extension=".txt";
	FILE*fOut;
	std::ofstream myfile;

	int cadhID, startNodeID, endNodeID, cellID, IDGap;

	cellID=floor(nCellRow/2);
	
	myfile.open ("zCadhSize.txt");
	 myfile << "Size: " << CadhPack.List.size() << "\n";
	 myfile.close();
	//std::cout << "Size:  " << CadhPack.List.size()  << std::endl;

	sprintf(name, "GAPS_POINTS_%s%s" , index_cluster, extension);				// 
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	fprintf(fOut, "STEP: %d ", t );
	fprintf (fOut, "\n");
	for (int i=0;i< ActiveGAP.size();++i)
	{
		IDGap=ActiveGAP[i];
		fprintf(fOut, "GAP: %d ", GAP.at(IDGap).IDGapType );
		fprintf (fOut, "\n");
		fprintf(fOut, "SHAPE:  %d ",GAP.at(IDGap).CurrentShapeIDList.size());
		for (int j=0; j<GAP.at(IDGap).CurrentShapeIDList.size(); j++)  
    		fprintf(fOut, "%d ", GAP.at(IDGap).CurrentShapeIDList[j]);	
    	
    	fprintf (fOut, "\n");
    	fprintf(fOut, "INIT POINT: %d  ; FINAL POINT: %d", GAP.at(IDGap).InitPoint, GAP.at(IDGap).LastPoint );
    	fprintf (fOut, "\n");
    	fprintf(fOut, "ID LIST: " );
    	fprintf(fOut, "%d :",GAP.at(IDGap).IDList.size());
    	for (int j=0; j<GAP.at(IDGap).IDList.size(); j++)  
    		fprintf(fOut, "%d ,", GAP.at(IDGap).IDList[j]);

    	fprintf (fOut, "\n");	
    }
    fprintf (fOut, "\n");
     fprintf (fOut, "\n");
    fclose(fOut);

  /*  
    sprintf(name, "GAPS_POINTS_INACTIVE_%s%s" , index_cluster, extension);				// 
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	fprintf(fOut, "STEP: %d ", t );
	fprintf (fOut, "\n");
	for (int i=0;i< InactiveGAP.size();++i)
	{
		fprintf(fOut, "GAP: %d ", i );
		fprintf(fOut, "%d ",InactiveGAP[i].CurrentShapeIDList.size());
		for (int j=0; j<InactiveGAP[i].CurrentShapeIDList.size(); j++)  
    		fprintf(fOut, "%d ", InactiveGAP[i].CurrentShapeIDList[j]);	
    	
    	fprintf (fOut, "\n");
    	fprintf(fOut, "INIT POINT: %d  ; FINAL POINT: %d", InactiveGAP[i].InitPoint[InactiveGAP[i].InitPoint.size()-1], InactiveGAP[i].LastPoint[InactiveGAP[i].LastPoint.size()-1] );
    	fprintf (fOut, "\n");
    	fprintf(fOut, "ID LIST: " );
    	fprintf(fOut, "%d :",InactiveGAP[i].IDList[InactiveGAP[i].IDList.size()-1].size());
    	for (int j=0; j<InactiveGAP[i].IDList[GAP[i].IDList.size()-1].size(); j++)  
    		fprintf(fOut, "%d ,", InactiveGAP[i].IDList[InactiveGAP[i].IDList.size()-1][j]);

    	fprintf (fOut, "\n");	
    }
    fprintf (fOut, "\n");
     fprintf (fOut, "\n");
    fclose(fOut);
*/

	sprintf(name, "MembraneDistance_%s%s" , index_cluster, extension);				// 
	fOut = fopen(name,"a");   //para reescribir se usa "a"
      

	
  if (nCellRow!=2)
  {
	 	//Side 0
	startNodeID=0;
	endNodeID=nAdhPointSide;

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		
		fprintf (fOut, "%4.1f ", MembDistance.at(i).distance);   	   	    		   	    		    
	}
 	   	    		   	    	
	fprintf(fOut, "\n");

	//Side 1
	startNodeID=endNodeID;
	endNodeID=endNodeID+nAdhPointSide;

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		
		fprintf (fOut, "%4.1f ", MembDistance.at(i).distance);   	   	    		   	    		    
	}
 	   	    		   	    	
	fprintf(fOut, "\n");

	

	//Side 2
	startNodeID=endNodeID;
	endNodeID=endNodeID+nAdhPointSide;
	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		
		fprintf (fOut, "%4.1f ", MembDistance.at(i).distance);   	   	    		   	    		    
	}
 	   	    		   	    	
	fprintf(fOut, "\n");

	

  }	

	//Side 3
  startNodeID=endNodeID;
	endNodeID=endNodeID+nAdhPointSide;
for (auto i=startNodeID; i< endNodeID; ++i )
	{
		
		fprintf (fOut, "%4.1f ", MembDistance.at(i).distance);   	   	    		   	    		    
	}
 	   	    		   	    	
	fprintf(fOut, "\n");

	

	//Side 4
	startNodeID=endNodeID;
	endNodeID=endNodeID+nAdhPointSide;
for (auto i=startNodeID; i< endNodeID; ++i )
	{
		
		fprintf (fOut, "%4.1f ", MembDistance.at(i).distance);   	   	    		   	    		    
	}
 	   	    		   	    	
	fprintf(fOut, "\n");

	
	//Side 5
	startNodeID=endNodeID;
	endNodeID=endNodeID+nAdhPointSide;
for (auto i=startNodeID; i< endNodeID; ++i )
	{
		
		fprintf (fOut, "%4.1f ", MembDistance.at(i).distance);   	   	    		   	    		    
	}
 	   	    		   	    	
	fprintf(fOut, "\n");

	
fclose(fOut);





///////////////////  I only use the  value of the last step
	sprintf(name, "GapTotalTime_%s%s" , index_cluster, extension);
	fOut = fopen(name,"a");   //para reescribir se usa "a"

  if (nCellRow!=2)
  {
	 	//Side 0
	startNodeID=GETVertexID(cellID, 0);
	endNodeID=GETVertexID(cellID, 1);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		//cadhID=Nodes.at(i)-> cadhPointer.at(0);
		fprintf (fOut, "%4.1f ", Nodes.at(i) -> gapTotTime);   	   	    		   	    		    
	}

	//cadhID=Nodes.at(endNodeID)-> cadhPointer.at(1);
//	fprintf (fOut, "%4.1f ", Nodes.at(endNodeID) -> gapTotTime);   	   	    		   	    	
	fprintf(fOut, "\n");

	//Side 1
	startNodeID=GETVertexID(cellID, 1);
	endNodeID=GETVertexID(cellID, 2);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
	//	cadhID=Nodes.at(i)-> cadhPointer.at(0);
		fprintf (fOut, "%4.1f ", Nodes.at(i) -> gapTotTime);   	   	    		   	    		    
	}

	//cadhID=Nodes.at(endNodeID)-> cadhPointer.at(0);
	//fprintf (fOut, "%4.1f ", Nodes.at(endNodeID) -> gapTotTime);   	   	    		   	    	
	fprintf(fOut, "\n");	

	

	//Side 2
	startNodeID=GETVertexID(cellID, 2);
	endNodeID=GETVertexID(cellID, 3);


	//cadhID=Nodes.at(startNodeID)-> cadhPointer.at(1);
	//fprintf (fOut, "%4.1f ", Nodes.at(startNodeID) -> gapTotTime);   	   	    		   	    	
	for (auto i=startNodeID; i< endNodeID; ++i )
	{
	//	cadhID=Nodes.at(i)-> cadhPointer.at(0);
		fprintf (fOut, "%4.1f ", Nodes.at(i) -> gapTotTime);   	   	    		   	    		    
	}

	//cadhID=Nodes.at(endNodeID)-> cadhPointer.at(1);
	//fprintf (fOut, "%4.1f ", Nodes.at(endNodeID) -> gapTotTime);   	   	    		   	    	
	fprintf(fOut, "\n");

  }	

	//Side 3
	startNodeID=GETVertexID(cellID, 3);
	endNodeID=GETVertexID(cellID, 4);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
	//	cadhID=Nodes.at(i)-> cadhPointer.at(0);
		fprintf (fOut, "%4.1f ", Nodes.at(i) -> gapTotTime);   	   	    		   	    		    
	}

	//cadhID=Nodes.at(endNodeID)-> cadhPointer.at(1);
	//fprintf (fOut, "%4.1f ", Nodes.at(endNodeID) -> gapTotTime);   	   	    		   	    	
	fprintf(fOut, "\n");


	//Side 4
	startNodeID=GETVertexID(cellID, 4);
	endNodeID=GETVertexID(cellID, 5);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
	//	cadhID=Nodes.at(i)-> cadhPointer.at(0);
		fprintf (fOut, "%4.1f ", Nodes.at(i) -> gapTotTime);   	   	    		   	    		    
	}

	//cadhID=Nodes.at(endNodeID)-> cadhPointer.at(0);
	//fprintf (fOut, "%4.1f ", Nodes.at(endNodeID) -> gapTotTime);   	   	    		   	    	
	fprintf(fOut, "\n");

	//Side 5
	startNodeID=GETVertexID(cellID, 5);
	endNodeID=startNodeID+nAdhPointSide;


	//cadhID=Nodes.at(startNodeID)-> cadhPointer.at(1);
	//fprintf (fOut, "%4.1f ", Nodes.at(startNodeID) -> gapTotTime);   	   	    		   	    	
	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		//cadhID=Nodes.at(i)-> cadhPointer.at(0);
		fprintf (fOut, "%4.1f ", Nodes.at(i) -> gapTotTime);   	   	    		   	    		    
	}
/*
	if (nCellRow!=2)
   {	
		endNodeID=GETVertexID(cellID, 0);
		//cadhID=Nodes.at(endNodeID)-> cadhPointer.at(1);
		fprintf (fOut, "%4.1f ", Nodes.at(endNodeID) -> gapTotTime);   
	}		   	    		   	    	
*/
	fprintf(fOut, "\n");
	fclose(fOut);


	sprintf(name, "InactiveAdhPoints_%s%s" , index_cluster, extension);   // Time to generate the first gap
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	fprintf (fOut, "Step: %d \n", t); 
	for (auto i=0; i< Pool.size(); i++ )
	{
		fprintf (fOut, "CELL: %d: ", i );   // Print mean traction and last traction	
		for (auto it=Pool[i].UnboundAdhPoint.begin(); it!=Pool[i].UnboundAdhPoint.end();)
		{
			
			fprintf (fOut, "%d ", *(it) );   // Print mean traction and last traction		    
			it++;
		}
		fprintf(fOut, "\n");
	}
	fclose(fOut);

	sprintf(name, "zCadherins_%s%s" , index_cluster, extension);   // Time to generate the first gap
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	fprintf(fOut, "\n");
		fprintf (fOut, "Step: %d \n", t); 
	
	for (auto it=0;it < CadhPack.List.size(); it++ )
	{

		fprintf (fOut, "Cadh: %d :", it);   // Print mean traction and last traction	
		fprintf (fOut, "%d :", CadhPack.List.at(it).firstAdhPoint);   // Print mean traction and last traction  
			fprintf (fOut, "%d  Active: %d \n:", CadhPack.List.at(it).secondAdhPoint, CadhPack.List.at(it).active);   // Print mean traction and last traction    
  
	}
    fclose(fOut);
	sprintf(name, "zCadherins_active_%s%s" , index_cluster, extension);   // Time to generate the first gap
	fOut = fopen(name,"a");   //para reescribir se usa "a"
		fprintf (fOut, "Step: %d \n", t); 
	
	for (auto it=0;it < CadhPack.ActiveCadh.size(); it++ )
	{

		fprintf (fOut, "Cadh: %d :", it);   // Print mean traction and last traction	
		fprintf (fOut, "%d :", CadhPack.List.at(CadhPack.ActiveCadh[it]).firstAdhPoint);   // Print mean traction and last traction  
			fprintf (fOut, "%d \n:", CadhPack.List.at(CadhPack.ActiveCadh[it]).secondAdhPoint);   // Print mean traction and last traction    
  
	}

 fclose(fOut);


 PrintGAPStep(index_cluster);
	
}


void Monolayer::printDataGlobal(char *index_cluster)
{

	char name[100];
	char orden[1024]; 
	char * extension=".txt";
	FILE*fOut;

/*
	sprintf(name, "InactiveAdhPoints_%s%s" , index_cluster, extension);   // Time to generate the first gap
	fOut = fopen(name,"a");   //para reescribir se usa "a"

	for (auto i=0; i< Pool.size(); i++ )
	{
		fprintf (fOut, "CELL: %d: ", i );   // Print mean traction and last traction	
		for(j=0; j < Pool[i].UnboundAdhPoint.size();j++ )
			fprintf (fOut, "%4.1f ", Pool[i].UnboundAdhPoint[j] );   // Print mean traction and last traction		    
		fprintf(fOut, "\n");
	}


	sprintf(name, "Cadherins_%s%s" , index_cluster, extension);   // Time to generate the first gap
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	int ct=0;
	for (auto it=CadhPack.List.begin(); it !=CadhPack.List.end(); )
	{

		fprintf (fOut, "Cadh: %d :", ct);   // Print mean traction and last traction	
		fprintf (fOut, "%d :", (*it).firstAdhPoint);   // Print mean traction and last traction  
				fprintf (fOut, "%d \n:", (*it).secondAdhPoint);   // Print mean traction and last traction    
  
		it++;
		ct++;
	}


	fclose(fOut);


	fclose(fOut);

	sprintf(name, "GapNumber_%s%s" , index_cluster, extension);			// Number of Gaps generated
	fOut = fopen(name,"a");   //para reescribir se usa "a"

	for (auto it=CadhPack.List.begin(); it !=CadhPack.List.end(); )
	{
		fprintf (fOut, "%d ", (*it).nGap);   // Print mean traction and last traction		    
		it++;
	}


	fclose(fOut);

	double aver;

	sprintf(name, "GapGenFreq_%s%s" , index_cluster, extension);   // averga of times that takes to generate the gap
	fOut = fopen(name,"a");   //para reescribir se usa "a"

	for (auto it=CadhPack.List.begin(); it !=CadhPack.List.end(); )
	{
		aver  = std::accumulate((*it).gapGenTime.begin(), (*it).gapGenTime.end(), double(0) ) / (*it).gapGenTime.size();
		fprintf (fOut, "%4.3f ", aver );   // Print mean traction and last traction		    
		it++;
	}


	fclose(fOut);


	sprintf(name, "GapGenFirst_%s%s" , index_cluster, extension);   // Time to generate the first gap
	fOut = fopen(name,"a");   //para reescribir se usa "a"

	for (auto it=CadhPack.List.begin(); it !=CadhPack.List.end(); )
	{
		fprintf (fOut, "%4.1f ", (*it).gapGenTime.at(0) );   // Print mean traction and last traction		    
		it++;
	}


	fclose(fOut);
*/
}

void Monolayer::PrintGAPStep(char *index_cluster)
{

	char name[100];
	char orden[1024]; 
	int ct=0;
	int IDGap;
	double vertexTime=double(0);
	double interfaceTime=double(0);
	double stressAverage;
	FILE*fOut;

	for (auto i=0; i< ActiveGAP.size(); ++i)
	{
		IDGap=ActiveGAP[i];
		stressAverage=GapStress(GAP.at(IDGap));
		if ( GAP.at(IDGap).vertex )
		{
			sprintf(name, "GapVertex_%s_%d.txt" , index_cluster, GAP.at(IDGap).ctGapsID);			// Number of Gaps generated
			fOut = fopen(name,"a");   //para reescribir se usa "a"
			ctVertexStep++;

		}
		else
		{
			sprintf(name, "GapInterface_%s_%d.txt" , index_cluster, GAP.at(IDGap).ctGapsID);			// Number of Gaps generated
			fOut = fopen(name,"a");   //para reescribir se usa "a"
			ctInterfaceStep++;
		}
		//for (ct=0; ct < GAP.at(IDGap).Area.size(); ++ct )
			//fprintf (fOut, "%d %d %d %d %4.4f %4.4f %4.4f %4.4f %4.4f\n ",step,  GAP.at(IDGap).nAdhPoint, GAP.at(IDGap).InitPoint, GAP.at(IDGap).LastPoint, 
			//			GAP.at(IDGap).Area, GAP.at(IDGap).Length, GAP.at(IDGap).HeightAver, GAP.at(IDGap).HeightMax, GAP.at(IDGap).DistVertex );   // Print mean traction and last traction		    
			fprintf (fOut, "%4.4f %4.4f %4.4f %d\n",step*TimeStep, GAP.at(IDGap).Area, stressAverage, GAP.at(IDGap).InitPoint);   // Print mean traction and last traction		    
		fclose(fOut);	
	}
}

double Monolayer::GapStress(GapType &GapAux )
{
	double forceTotal=double(0), forcemod, lCurrent, lBalance;
	int pointA, pointB,startAux, pointAux;
	startAux=GapAux.InitPoint;
	CellIDfromPoint(startAux);
	pointAux=GETVertexID(CellIDfromPoint(startAux), 0);

	for (int i=0;  i< GapAux.IDList.size()+1; ++i)
	{
		pointA=NormalizedID(startAux+i,pointAux); 
		pointB=NormalizedID(pointA-1, pointAux); 

		lCurrent=NodesDist(pointA, pointB);
		lBalance=sepAdhPoint;
		//EACurrent=EABalance*lBalance/lCurrent;
		// Using EA as K
		forcemod=membEA*(lCurrent-lBalance);
		forceTotal+=forcemod;
	}
	forceTotal=forceTotal/double((GapAux.IDList.size()+1) );
	return forceTotal;

}


void Monolayer::PrintGlobalStat(char *index_cluster)
{
	char name[100];
	char orden[1024]; 
	
	FILE*fOut;

	sprintf(name, "IndexFile.txt" , index_cluster);			// Number of Gaps generated
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	fprintf (fOut, "%s", index_cluster);
	fclose(fOut);

	sprintf(name, "GapStat_%s.txt" , index_cluster);			// Number of Gaps generated
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	// Total vertex time_steps, interface, number of vertex, number of interface, steos 
	//fprintf (fOut, "%d %d %d %d %d %d\n ",ctVertexStep,  ctInterfaceStep,  ctVertexPrint, ctInterfacePrint, ctVertexFirst, ctInterfaceFirst );
	fprintf (fOut, "%4.6f %4.6f %4.6f %4.6f %4.6f %4.6f\n ",ctVertexStep*TimeStep,  ctInterfaceStep*TimeStep,  
					ctVertexPrint/(TimeStep*TimeSlots), ctInterfacePrint/(TimeStep*TimeSlots), 
					ctVertexFirst*TimeStep, ctInterfaceFirst*TimeStep );
	fclose(fOut);	

	if (generalRupture)
	{	
		sprintf(name, "BIG_GAP_%s.txt" , index_cluster);			// Number of Gaps generated
		fOut = fopen(name,"a");   //para reescribir se usa "a"
		fclose(fOut);
	}
	
}

void Monolayer::PrintStepStat(char *index_cluster)
{
	char name[100];
	char orden[1024]; 
	int IDAux;
	
	FILE*fOut;
	int ctCadh=0;


	for (int j=0; j < CadhPack.ActiveCadh.size(); ++j)  
	{
		IDAux=CadhPack.ActiveCadh.at(j);
		ctCadh=ctCadh+CadhPack.List.at(IDAux).cadherinNumber;

	}


	sprintf(name, "Global_Adhesion_%s.txt" , index_cluster);			// Number of Gaps generated
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	// Total vertex time_steps, interface, number of vertex, number of interface, steos 
	//fprintf (fOut, "%d %d %d %d %d %d\n ",ctVertexStep,  ctInterfaceStep,  ctVertexPrint, ctInterfacePrint, ctVertexFirst, ctInterfaceFirst );
	fprintf (fOut, "%d %d\n", ctCadh, CadhPack.ActiveCadh.size() );
	fclose(fOut);	
	
}



void Monolayer::printDataGlobalVertex(char *index_cluster)
{

	char name[100];
	char orden[1024]; 
	char *extension=".txt";
	FILE*fOut;


	int cellID, startNodeID, endNodeID, cadhID;


	cellID=floor(nCellRow/2);



	

	sprintf(name, "GapNumber_%s%s" , index_cluster, extension);			// Number of Gaps generated
	fOut = fopen(name,"a");   //para reescribir se usa "a"

 	//Side 0


  if (nCellRow!=2)
  {	
	startNodeID=GETVertexID(cellID, 0);
	endNodeID=GETVertexID(cellID, 1);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
	//	cadhID=Nodes.at(i)-> cadhPointer.at(0);
		fprintf (fOut, "%d ", Nodes.at(i) -> nGap);   // Print mean traction and last traction		    
	}

	//cadhID=Nodes.at(endNodeID)-> cadhPointer.at(1);
	//fprintf (fOut, "%d ", Nodes.at(endNodeID) -> nGap);   // Print mean traction and last traction	
	fprintf(fOut, "\n");

	//Side 1
	startNodeID=GETVertexID(cellID, 1);
	endNodeID=GETVertexID(cellID, 2);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
	//	cadhID=Nodes.at(i)-> cadhPointer.at(0);
		fprintf (fOut, "%d ", Nodes.at(i) -> nGap);   // Print mean traction and last traction		    
	}

	//cadhID=Nodes.at(endNodeID)-> cadhPointer.at(0);
	//fprintf (fOut, "%d ", Nodes.at(endNodeID) -> nGap);   // Print mean traction and last traction	
	fprintf(fOut, "\n");	

	

	//Side 2
	startNodeID=GETVertexID(cellID, 2);
	endNodeID=GETVertexID(cellID, 3);


	// cadhID=Nodes.at(startNodeID)-> cadhPointer.at(1);
	//fprintf (fOut, "%d ", Nodes.at(startNodeID) -> nGap);   // Print mean traction and last traction	
	for (auto i=startNodeID; i< endNodeID; ++i )
	{
	//	// cadhID=Nodes.at(i)-> cadhPointer.at(0);
		fprintf (fOut, "%d ", Nodes.at(i) -> nGap);   // Print mean traction and last traction		    
	}

	// cadhID=Nodes.at(endNodeID)-> cadhPointer.at(1);
	//fprintf (fOut, "%d ", Nodes.at(endNodeID) -> nGap);   // Print mean traction and last traction	
	fprintf(fOut, "\n");
  }	

	//Side 3
	startNodeID=GETVertexID(cellID, 3);
	endNodeID=GETVertexID(cellID, 4);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		// cadhID=Nodes.at(i)-> cadhPointer.at(0);
		fprintf (fOut, "%d ", Nodes.at(i) -> nGap);   // Print mean traction and last traction		    
	}

	// cadhID=Nodes.at(endNodeID)-> cadhPointer.at(1);
	//fprintf (fOut, "%d ", Nodes.at(endNodeID) -> nGap);   // Print mean traction and last traction	
	fprintf(fOut, "\n");


	//Side 4
	startNodeID=GETVertexID(cellID, 4);
	endNodeID=GETVertexID(cellID, 5);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		// cadhID=Nodes.at(i)-> cadhPointer.at(0);
		fprintf (fOut, "%d ", Nodes.at(i) -> nGap);   // Print mean traction and last traction		    
	}

	// cadhID=Nodes.at(endNodeID)-> cadhPointer.at(0);
	//fprintf (fOut, "%d ", Nodes.at(endNodeID) -> nGap);   // Print mean traction and last traction	
	fprintf(fOut, "\n");

	//Side 5
	startNodeID=GETVertexID(cellID, 5);
	endNodeID=startNodeID+nAdhPointSide;


	// cadhID=Nodes.at(startNodeID)-> cadhPointer.at(1);
//	fprintf (fOut, "%d ", Nodes.at(startNodeID) -> nGap);   // Print mean traction and last traction	
	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		// cadhID=Nodes.at(i)-> cadhPointer.at(0);
		fprintf (fOut, "%d ", Nodes.at(i) -> nGap);   // Print mean traction and last traction		    
	}

 
	fclose(fOut);	



	//////////////////////////////

	sprintf(name, "GapGenFreq_%s%s" , index_cluster, extension);			// Number of Gaps generated
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	double aver;
 	//Side 0



  if (nCellRow!=2)
  {	
	startNodeID=GETVertexID(cellID, 0);
	endNodeID=GETVertexID(cellID, 1);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		// cadhID=Nodes.at(i)-> cadhPointer.at(0);
		aver  = std::accumulate(Nodes.at(i) -> gapGenTime.begin(), 
								Nodes.at(i) -> gapGenTime.end(), double(0) ) / Nodes.at(i) -> gapGenTime.size();
		fprintf (fOut, "%4.3f ", aver );   // Print mean traction and last traction		    
		    
	}

	// cadhID=Nodes.at(endNodeID)-> cadhPointer.at(1);
	//aver  = std::accumulate(Nodes.at(endNodeID) -> gapGenTime.begin(), 
	//							Nodes.at(endNodeID) -> gapGenTime.end(), double(0) ) / Nodes.at(endNodeID) -> gapGenTime.size();
	//fprintf (fOut, "%4.3f ", aver );   // Print mean traction and last traction		    
	fprintf(fOut, "\n");

	//Side 1
	startNodeID=GETVertexID(cellID, 1);
	endNodeID=GETVertexID(cellID, 2);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		// cadhID=Nodes.at(i)-> cadhPointer.at(0);
		aver  = std::accumulate(Nodes.at(i) -> gapGenTime.begin(), 
								Nodes.at(i) -> gapGenTime.end(), double(0) ) / Nodes.at(i) -> gapGenTime.size();
		fprintf (fOut, "%4.3f ", aver );   // Print mean traction and last traction		    
	    
	}

	// cadhID=Nodes.at(endNodeID)-> cadhPointer.at(0);
	aver  = std::accumulate(Nodes.at(endNodeID) -> gapGenTime.begin(), 
								Nodes.at(endNodeID) -> gapGenTime.end(), double(0) ) / Nodes.at(endNodeID) -> gapGenTime.size();
	//fprintf (fOut, "%4.3f ", aver );   // Print mean traction and last traction	
	fprintf(fOut, "\n");	

	

	//Side 2
	startNodeID=GETVertexID(cellID, 2);
	endNodeID=GETVertexID(cellID, 3);


	// cadhID=Nodes.at(startNodeID)-> cadhPointer.at(1);
	aver  = std::accumulate(Nodes.at(startNodeID) -> gapGenTime.begin(), 
								Nodes.at(startNodeID) -> gapGenTime.end(), double(0) ) / Nodes.at(startNodeID) -> gapGenTime.size();
//	fprintf (fOut, "%4.3f ", aver );   // Print mean traction and last traction	
	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		// cadhID=Nodes.at(i)-> cadhPointer.at(0);
		aver  = std::accumulate(Nodes.at(i) -> gapGenTime.begin(), 
								Nodes.at(i) -> gapGenTime.end(), double(0) ) / Nodes.at(i) -> gapGenTime.size();
		fprintf (fOut, "%4.3f ", aver );   // Print mean traction and last traction	
	}

	// cadhID=Nodes.at(endNodeID)-> cadhPointer.at(1);
	aver  = std::accumulate(Nodes.at(endNodeID) -> gapGenTime.begin(), 
								Nodes.at(endNodeID) -> gapGenTime.end(), double(0) ) / Nodes.at(endNodeID) -> gapGenTime.size();
//	fprintf (fOut, "%4.3f ", aver );   // Print mean traction and last traction	
	fprintf(fOut, "\n");
  }
	//Side 3
	startNodeID=GETVertexID(cellID, 3);
	endNodeID=GETVertexID(cellID, 4);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		// cadhID=Nodes.at(i)-> cadhPointer.at(0);
		aver  = std::accumulate(Nodes.at(i) -> gapGenTime.begin(), 
								Nodes.at(i) -> gapGenTime.end(), double(0) ) / Nodes.at(i) -> gapGenTime.size();
		fprintf (fOut, "%4.3f ", aver );   // Print mean traction and last traction	
	}

	// cadhID=Nodes.at(endNodeID)-> cadhPointer.at(1);
	aver  = std::accumulate(Nodes.at(endNodeID) -> gapGenTime.begin(), 
								Nodes.at(endNodeID) -> gapGenTime.end(), double(0) ) / Nodes.at(endNodeID) -> gapGenTime.size();
	//fprintf (fOut, "%4.3f ", aver );   // Print mean traction and last traction	
	fprintf(fOut, "\n");


	//Side 4
	startNodeID=GETVertexID(cellID, 4);
	endNodeID=GETVertexID(cellID, 5);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		// cadhID=Nodes.at(i)-> cadhPointer.at(0);
		aver  = std::accumulate(Nodes.at(i) -> gapGenTime.begin(), 
								Nodes.at(i) -> gapGenTime.end(), double(0) ) / Nodes.at(i) -> gapGenTime.size();
		fprintf (fOut, "%4.3f ", aver );   // Print mean traction and last traction	
	}

	// cadhID=Nodes.at(endNodeID)-> cadhPointer.at(0);
	aver  = std::accumulate(Nodes.at(endNodeID) -> gapGenTime.begin(), 
								Nodes.at(endNodeID) -> gapGenTime.end(), double(0) ) / Nodes.at(endNodeID) -> gapGenTime.size();
//fprintf (fOut, "%4.3f ", aver );   // Print mean traction and last traction	
	fprintf(fOut, "\n");

	//Side 5
	startNodeID=GETVertexID(cellID, 5);
	endNodeID=startNodeID+nAdhPointSide;


	// cadhID=Nodes.at(startNodeID)-> cadhPointer.at(1);
	aver  = std::accumulate(Nodes.at(startNodeID) -> gapGenTime.begin(), 
								Nodes.at(startNodeID) -> gapGenTime.end(), double(0) ) / Nodes.at(startNodeID) -> gapGenTime.size();
	//fprintf (fOut, "%4.3f ", aver );   // Print mean traction and last traction	
	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		// cadhID=Nodes.at(i)-> cadhPointer.at(0);
		aver  = std::accumulate(Nodes.at(i) -> gapGenTime.begin(), 
								Nodes.at(i) -> gapGenTime.end(), double(0) ) / Nodes.at(i) -> gapGenTime.size();
		fprintf (fOut, "%4.3f ", aver );   // Print mean traction and last traction	
	}

/*
  if (nCellRow!=2)
  {
	endNodeID=GETVertexID(cellID, 0);
	// cadhID=Nodes.at(endNodeID)-> CadhPointer.at(1);
	aver  = std::accumulate(Nodes.at(endNodeID) -> gapGenTime.begin(), 
								Nodes.at(endNodeID) -> gapGenTime.end(), double(0) ) / Nodes.at(endNodeID) -> gapGenTime.size();
	fprintf (fOut, "%4.3f ", aver );   // Print mean traction and last traction	
//	fprintf (fOut, "%d ", Nodes.at(i) -> nGap);   // Print mean traction and last traction	
  }
  */
	fclose(fOut);	



	////////////


	sprintf(name, "GapGenFirst_%s%s" , index_cluster, extension);			// Number of Gaps generated
	fOut = fopen(name,"a");   //para reescribir se usa "a"

 	//Side 0

  if (nCellRow!=2)
  {	
	startNodeID=GETVertexID(cellID, 0);
	endNodeID=GETVertexID(cellID, 1);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		// cadhID=Nodes.at(i)-> CadhPointer.at(0);
		fprintf (fOut, "%4.1f ", Nodes.at(i) -> gapGenTime.at(0) );   // Print mean traction and last traction		    		    
	}

	// cadhID=Nodes.at(endNodeID)-> CadhPointer.at(1);
	//fprintf (fOut, "%4.1f ", Nodes.at(endNodeID) -> gapGenTime.at(0) );   // Print mean traction and last traction		    	
	fprintf(fOut, "\n");

	//Side 1
	startNodeID=GETVertexID(cellID, 1);
	endNodeID=GETVertexID(cellID, 2);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		// cadhID=Nodes.at(i)-> CadhPointer.at(0);
		fprintf (fOut, "%4.1f ", Nodes.at(i) -> gapGenTime.at(0) );   // Print mean traction and last traction		    		    
	}

	// cadhID=Nodes.at(endNodeID)-> CadhPointer.at(0);
	//fprintf (fOut, "%4.1f ", Nodes.at(endNodeID) -> gapGenTime.at(0) );   // Print mean traction and last traction		    	
	fprintf(fOut, "\n");	

	

	//Side 2
	startNodeID=GETVertexID(cellID, 2);
	endNodeID=GETVertexID(cellID, 3);


	// cadhID=Nodes.at(startNodeID)-> CadhPointer.at(1);
	//fprintf (fOut, "%4.1f ", Nodes.at(startNodeID) -> gapGenTime.at(0) );   // Print mean traction and last traction		    	
	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		// cadhID=Nodes.at(i)-> CadhPointer.at(0);
		fprintf (fOut, "%4.1f ", Nodes.at(i) -> gapGenTime.at(0) );   // Print mean traction and last traction		    		    
	}

	// cadhID=Nodes.at(endNodeID)-> CadhPointer.at(1);
	//fprintf (fOut, "%4.1f ", Nodes.at(endNodeID) -> gapGenTime.at(0) );   // Print mean traction and last traction		    	
	fprintf(fOut, "\n");
   }
	//Side 3
	startNodeID=GETVertexID(cellID, 3);
	endNodeID=GETVertexID(cellID, 4);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		// cadhID=Nodes.at(i)-> CadhPointer.at(0);
		fprintf (fOut, "%4.1f ", Nodes.at(i) -> gapGenTime.at(0) );   // Print mean traction and last traction		    		    
	}

	// cadhID=Nodes.at(endNodeID)-> CadhPointer.at(1);
	//fprintf (fOut, "%4.1f ", Nodes.at(endNodeID) -> gapGenTime.at(0) );   // Print mean traction and last traction		    	
	fprintf(fOut, "\n");


	//Side 4
	startNodeID=GETVertexID(cellID, 4);
	endNodeID=GETVertexID(cellID, 5);

	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		// cadhID=Nodes.at(i)-> CadhPointer.at(0);
		fprintf (fOut, "%4.1f ", Nodes.at(i) -> gapGenTime.at(0) );   // Print mean traction and last traction		    		    
	}

	// cadhID=Nodes.at(endNodeID)-> CadhPointer.at(0);
	//fprintf (fOut, "%4.1f ", Nodes.at(endNodeID) -> gapGenTime.at(0) );   // Print mean traction and last traction		    	
	fprintf(fOut, "\n");

	//Side 5
	startNodeID=GETVertexID(cellID, 5);
	endNodeID=startNodeID+nAdhPointSide;


	// cadhID=Nodes.at(startNodeID)-> CadhPointer.at(1);
	//fprintf (fOut, "%4.1f ", Nodes.at(startNodeID) -> gapGenTime.at(0) );   // Print mean traction and last traction		    	
	for (auto i=startNodeID; i< endNodeID; ++i )
	{
		//// cadhID=Nodes.at(i)-> CadhPointer.at(0);
		fprintf (fOut, "%4.1f ", Nodes.at(i) -> gapGenTime.at(0) );   // Print mean traction and last traction		    		    
	}
/*
  if (nCellRow!=2)
  {
	endNodeID=GETVertexID(cellID, 0);
	// cadhID=Nodes.at(endNodeID)-> CadhPointer.at(1);
	fprintf (fOut, "%4.1f ", Nodes.at(endNodeID) -> gapGenTime.at(0) );   // Print mean traction and last traction		    	
  }*/
	fclose(fOut);	

}

void Monolayer::printVTK(int ctParaview)
  {
	int j,c, nConnectivity;
	int IDAux;
	char name[50];
	char pref[]="Cell_";	
	char extension[]=".vtk";
	std::map<int, std::shared_ptr<MemAdhPoint>>::iterator it;

	Eigen::VectorXd FInt;
	FInt=Eigen::VectorXd::Zero(2*meshSize);
	CalculateForceIntern(FInt);
	ApplyBdryConditions(FInt);

	std::vector<double> vDir;
	vDir.resize(2);

	nConnectivity=ConnectivityMemb.at(0).size()+ConnectivityBdry.at(0).size();
		
	FILE*fOut;
	//Imprime archivo vtk
  	sprintf(name, "%s%d%s", pref, ctParaview , extension);
 	fOut = fopen(name,"w");   //para reescribir se usa "a"
	fprintf (fOut, "# vtk DataFile Version 4.1\n");
	fprintf (fOut, "Cell data\n");
	fprintf (fOut, "ASCII\n");
	fprintf (fOut, "DATASET POLYDATA\n");
	fprintf (fOut, "POINTS %d double\n",Nodes.size());
	
	for (it=Nodes.begin(); it !=Nodes.end(); ++it) 	fprintf (fOut, "%5.4f %5.4f %5.4f\n", 
													(it->second)->Coord.at(0), (it->second)->Coord.at(1), double(0) );

	fprintf (fOut, "\n");
	fprintf (fOut, "LINES %d %d\n",nConnectivity, 3*nConnectivity);
  	for (c = 0; c < ConnectivityMemb.at(0).size(); c++)  
    fprintf(fOut, "%d %d %d \n", 2, ConnectivityMemb.at(0).at(c),ConnectivityMemb.at(1).at(c));
	for (c = 0; c < ConnectivityBdry.at(0).size(); c++) 
    fprintf(fOut, "%d %d %d \n", 2, ConnectivityBdry.at(0).at(c),ConnectivityBdry.at(1).at(c));
		
	fprintf (fOut, "\n");	
	fprintf (fOut, "POINT_DATA %d\n",Nodes.size());
	fprintf (fOut, "SCALARS Bound double\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	
	for (it=Nodes.begin(); it !=Nodes.end(); ++it)  fprintf (fOut, "%d ", 	(it->second)->get_bound() );
	fprintf (fOut, "\n");

	fprintf (fOut, "SCALARS ID int\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	c=0;
	for (it=Nodes.begin(); it !=Nodes.end(); ++it)  
	{
		fprintf (fOut, "%d ", c);
		c++;
	}
	fprintf (fOut, "\n");
	fprintf (fOut, "SCALARS GAP_USED double\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	
	for (it=Nodes.begin(); it !=Nodes.end(); ++it)  fprintf (fOut, "%d ", 	(it->second)->GETGapIncluded() );
	fprintf (fOut, "\n");

	fprintf (fOut, "\n");
	fprintf (fOut, "VECTORS ForceBrowninan(pN) double\n");
	
	for (c=0; c< meshSize;++c )  
		fprintf (fOut, "%4.6f %4.6f 0.0\n", ForceStepBr(2*c)*1000.,ForceStepBr(2*c+1)*1000. );

	fprintf (fOut, "\n");



	fprintf (fOut, "VECTORS ForceInt(pN) double\n");
	
	for (c=0; c< meshSize;++c )  
		fprintf (fOut, "%4.6f %4.6f 0.0\n", FInt(2*c)*1000.,FInt(2*c+1)*1000. );

	fprintf (fOut, "\n");

	FInt=Eigen::VectorXd::Zero(2*meshSize);
	CalculateForceBdry(FInt);

	fprintf (fOut, "VECTORS ForceRadius(pN) double\n");
	
	for (c=0; c< meshSize;++c )  
		fprintf (fOut, "%4.6f %4.6f 0.0\n", FInt(2*c)*1000.,FInt(2*c+1)*1000. );

	fprintf (fOut, "\n");

	FInt=Eigen::VectorXd::Zero(2*meshSize);
	CalculateForceMemb(FInt);

	fprintf (fOut, "VECTORS ForceMemb(pN) double\n");
	
	for (c=0; c< meshSize;++c )  
		fprintf (fOut, "%4.6f %4.6f 0.0\n", FInt(2*c)*1000.,FInt(2*c+1)*1000. );

	fprintf (fOut, "\n");





	FInt=Eigen::VectorXd::Zero(2*meshSize);
//	RepulsionCellInCenter(FInt);
	RepulsionMembraneGlobalForce(FInt);

	fprintf (fOut, "VECTORS Repulsion(pN) double\n");
	
	for (c=0; c< meshSize;++c )  
		fprintf (fOut, "%4.6f %4.6f 0.0\n", FInt(2*c)*1000.,FInt(2*c+1)*1000. );

	fprintf (fOut, "\n");



	FInt=Eigen::VectorXd::Zero(2*meshSize);
	CalculateForceBendingMemb(FInt);

	fprintf (fOut, "VECTORS ForceBendingMemb(pN) double\n");
	
	for (c=0; c< meshSize;++c )  
		fprintf (fOut, "%4.6f %4.6f 0.0\n", FInt(2*c)*1000.,FInt(2*c+1)*1000. );

	fprintf (fOut, "\n");

	FInt=Eigen::VectorXd::Zero(2*meshSize);
	CalculateForceCadh(FInt);

	fprintf (fOut, "VECTORS ForceCadh(pN) double\n");
	
	for (c=0; c< meshSize;++c )  
		fprintf (fOut, "%4.6f %4.6f 0.0\n", FInt(2*c)*1000.,FInt(2*c+1)*1000. );

	fprintf (fOut, "\n");

	FInt=Eigen::VectorXd::Zero(2*meshSize);
	CalculateForceBendingCadh(FInt);

	fprintf (fOut, "VECTORS ForceBendingCadh(pN) double\n");
	
	for (c=0; c< meshSize;++c )  
		fprintf (fOut, "%4.6f %4.6f 0.0\n", FInt(2*c)*1000.,FInt(2*c+1)*1000. );

	fprintf (fOut, "\n");


	FInt=Eigen::VectorXd::Zero(2*meshSize);
	CalculateAngleInit(FInt);

	fprintf (fOut, "VECTORS AngleInit(pN) double\n");
	
	for (c=0; c< meshSize;++c )  
		fprintf (fOut, "%4.6f %4.6f 0.0\n", FInt(2*c),FInt(2*c+1) );

	fprintf (fOut, "\n");

	fprintf (fOut, "VECTORS ForceStep(pN) double\n");
	
	for (c=0; c< meshSize;++c )  
		fprintf (fOut, "%4.6f %4.6f 0.0\n", ForceStep(2*c)*1000.,ForceStep(2*c+1)*1000. );

	fprintf (fOut, "\n");

	
	// FInt=Eigen::VectorXd::Zero(2*meshSize);
	// SFOverlapForceVTK(FInt);

	// fprintf (fOut, "VECTORS SFOverlapForce(pN) double\n");
	
	// for (c=0; c< meshSize;++c )  
	// 	fprintf (fOut, "%4.6f %4.6f 0.0\n", FInt(2*c)*1000.,FInt(2*c+1)*1000. );

	// fprintf (fOut, "\n");


	

	
	fprintf (fOut, "VECTORS DirectionsImpeded float\n");
	
	for (c=0; c< meshSize;++c )  
		fprintf (fOut, "%d.0 %d.0 0.0\n", BdryPrint[2*c] ,BdryPrint[2*c+1] );

	fprintf (fOut, "\n");

	fprintf (fOut, "CELL_DATA %d\n",nConnectivity);

	fprintf (fOut, "scalars Stress double\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	
	for (c = 0; c < ConnectivityMemb.at(0).size(); c++) 
    	fprintf (fOut, "%4.4f ", CalculateStress(FInt, membEA, NodesDistInit( ConnectivityMemb.at(0).at(c),ConnectivityMemb.at(1).at(c) ), 
							ConnectivityMemb.at(0).at(c), ConnectivityMemb.at(1).at(c)) );
	for (c = 0; c < ConnectivityBdry.at(0).size(); c++) 
		fprintf (fOut, "%4.4f ", CalculateStress(FInt, bdryEA, SFlength[ConnectivityBdry.at(1).at(c)],
						ConnectivityBdry.at(0).at(c), ConnectivityBdry.at(1).at(c))  );
	
	fprintf (fOut, "\n");


	fprintf (fOut, "scalars Balance_Length_Inc double\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	
	for (c = 0; c < ConnectivityMemb.at(0).size(); c++) 
    	fprintf (fOut, "%4.4f ", double(0));
	for (c = 0; c < ConnectivityBdry.at(0).size(); c++) 
		fprintf (fOut, "%4.4f ", SFlength[ConnectivityBdry.at(1).at(c)]-SFlengthInit[ConnectivityBdry.at(1).at(c)]	);
	
	fprintf (fOut, "\n");


	fprintf (fOut, "scalars Balance_Length_Current double\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	
	for (c = 0; c < ConnectivityMemb.at(0).size(); c++) 
    	fprintf (fOut, "%4.4f ", double(aristLength));
	for (c = 0; c < ConnectivityBdry.at(0).size(); c++) 
		fprintf (fOut, "%4.4f ", SFlength[ConnectivityBdry.at(1).at(c)]	);
	
	fprintf (fOut, "\n");


	

	fprintf (fOut, "scalars Angle double\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	
	for (c = 0; c < ConnectivityMemb.at(0).size(); c++) 
    	fprintf (fOut, "%4.2f ", NodesAngle(ConnectivityMemb.at(0).at(c),ConnectivityMemb.at(1).at(c))*180/PI);
	for (c = 0; c < ConnectivityBdry.at(0).size(); c++) 
		fprintf (fOut, "%4.2f ", NodesAngle(ConnectivityBdry.at(0).at(c),ConnectivityBdry.at(1).at(c))*180/PI);
	
	fprintf (fOut, "\n");


	
	

	fprintf (fOut, "VECTORS AngleInit double\n");
	//fprintf (fOut, "LOOKUP_TABLE default\n");
	
	for (c = 0; c < ConnectivityMemb.at(0).size(); c++) 
	{
		
    	//fprintf (fOut, "%4.2f ", NodesAngleInit(ConnectivityMemb.at(0).at(c),ConnectivityMemb.at(1).at(c))*180/PI);
  	    fprintf (fOut, "%4.2f ", AngleBalanceCenter(ConnectivityMemb.at(0).at(c),ConnectivityMemb.at(1).at(c))*180/PI);
	}
	for (c = 0; c < ConnectivityBdry.at(0).size(); c++) 
	{	//fprintf (fOut, "%4.2f ", NodesAngleInit(ConnectivityBdry.at(0).at(c),ConnectivityBdry.at(1).at(c))*180/PI);
	  	fprintf (fOut, "%4.2f ", AngleBalanceCenter(ConnectivityMemb.at(0).at(c),ConnectivityMemb.at(1).at(c))*180/PI);
	
	}
	fprintf (fOut, "\n");




	fclose(fOut);
	int a,b;

	char pref5[]="ClosestP_";	
		
	//Imprime archivo vtk
  	sprintf(name, "%s%d%s", pref5, ctParaview , extension);
 	fOut = fopen(name,"w");   //para reescribir se usa "a"
	fprintf (fOut, "# vtk DataFile Version 4.1\n");
	fprintf (fOut, "Cadherin data\n");
	fprintf (fOut, "ASCII\n");
	fprintf (fOut, "DATASET POLYDATA\n");
	fprintf (fOut, "POINTS %d double\n",MembDistance.size()*2);

	
	
	for (int i=0; i < MembDistance.size(); i++)	
	{
		a=MembDistance[i].refID;
		b=MembDistance[i].pointMinID;
		//fprintf (fOut, "%5.4f %5.4f %5.4f\n",(Nodes.find(*it2)-> second)->Coord.at(0), (Nodes.find(*it2)-> second)->Coord.at(1), double(0) );
		fprintf (fOut, "%5.4f %5.4f %5.4f\n",(Nodes.at( a ))->Coord.at(0), (Nodes.at( a))->Coord.at(1), double(0) );
		fprintf (fOut, "%5.4f %5.4f %5.4f\n",(Nodes.at( b ))->Coord.at(0), (Nodes.at( b ))->Coord.at(1), double(0) );
	}

	fprintf (fOut, "\n");
	fprintf (fOut, "LINES %d %d\n",MembDistance.size(), 3*MembDistance.size());
  	for (c = 0; c < MembDistance.size(); c++) 
    fprintf(fOut, "%d %d %d \n", 2, 2*c,2*c+1);


	fprintf (fOut, "\n");

	fprintf (fOut, "CELL_DATA %d\n",MembDistance.size());

	fprintf (fOut, "SCALARS length double\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	
	 for (int i = 0; i < MembDistance.size(); i++)  
	 {
	 		
	 	a=MembDistance[i].refID;
		b=MembDistance[i].pointMinID;
	 	fprintf (fOut, "%4.4f ", 	NodesDist(a,b) );
	 }

	fprintf (fOut, "\n");
	fprintf (fOut, "\n");

	fprintf (fOut, "VECTORS ForceFinal(pN) double\n");
	
		
	for (int i=0; i < MembDistance.size(); i++)	
	{
		a=MembDistance[i].refID;
		b=MembDistance[i].pointMinID; 
		fprintf (fOut, "%d %d 0.0\n", a,b );
	}

	fprintf (fOut, "\n");





		
	fprintf (fOut, "\n");	
//	fprintf (fOut, "POINT_DATA %d\n",CadhPack.ActiveCadh.size()*2);
//	fprintf (fOut, "scalars Bound int\n");
//	fprintf (fOut, "LOOKUP_TABLE default\n");


	fclose(fOut);


	char pref7[]="ClosestPGlob_";	
		
	//Imprime archivo vtk
  	sprintf(name, "%s%d%s", pref7, ctParaview , extension);
 	fOut = fopen(name,"w");   //para reescribir se usa "a"
	fprintf (fOut, "# vtk DataFile Version 4.1\n");
	fprintf (fOut, "Cadherin data\n");
	fprintf (fOut, "ASCII\n");
	fprintf (fOut, "DATASET POLYDATA\n");
	fprintf (fOut, "POINTS %d double\n",MembDistanceGlobal.size()*2);

	
	
	for (int i=0; i < MembDistanceGlobal.size(); i++)	
	{
		a=MembDistanceGlobal[i].firstPoint;
		b=MembDistanceGlobal[i].secondPoint;
		//fprintf (fOut, "%5.4f %5.4f %5.4f\n",(Nodes.find(*it2)-> second)->Coord.at(0), (Nodes.find(*it2)-> second)->Coord.at(1), double(0) );
		fprintf (fOut, "%5.4f %5.4f %5.4f\n",(Nodes.at( a ))->Coord.at(0), (Nodes.at( a))->Coord.at(1), double(0) );
		fprintf (fOut, "%5.4f %5.4f %5.4f\n",(Nodes.at( b ))->Coord.at(0), (Nodes.at( b ))->Coord.at(1), double(0) );
	}

	fprintf (fOut, "\n");
	fprintf (fOut, "LINES %d %d\n",MembDistanceGlobal.size(), 3*MembDistanceGlobal.size());
  	for (c = 0; c < MembDistanceGlobal.size(); c++) 
    fprintf(fOut, "%d %d %d \n", 2, 2*c,2*c+1);


	fprintf (fOut, "\n");

	fprintf (fOut, "CELL_DATA %d\n",MembDistanceGlobal.size());

	fprintf (fOut, "SCALARS length double\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	
	 for (int i = 0; i < MembDistanceGlobal.size(); i++)  
	 {
	 		
	 	a=MembDistanceGlobal[i].firstPoint;
		b=MembDistanceGlobal[i].secondPoint;
	 	fprintf (fOut, "%4.4f ", 	NodesDist(a,b) );
	 }

	fprintf (fOut, "\n");


	




		
	fprintf (fOut, "\n");	
//	fprintf (fOut, "POINT_DATA %d\n",CadhPack.ActiveCadh.size()*2);
//	fprintf (fOut, "scalars Bound int\n");
//	fprintf (fOut, "LOOKUP_TABLE default\n");


	fclose(fOut);

	char pref2[]="CadherinBound_";	
		
	//Imprime archivo vtk
  	sprintf(name, "%s%d%s", pref2, ctParaview , extension);
 	fOut = fopen(name,"w");   //para reescribir se usa "a"
	fprintf (fOut, "# vtk DataFile Version 4.1\n");
	fprintf (fOut, "Cadherin data\n");
	fprintf (fOut, "ASCII\n");
	fprintf (fOut, "DATASET POLYDATA\n");
	fprintf (fOut, "POINTS %d double\n",CadhPack.ActiveCadh.size()*2);
	
	for (auto it2=CadhPack.ActiveCadh.begin(); it2 !=CadhPack.ActiveCadh.end(); ++it2) 	
	{
		//fprintf (fOut, "%5.4f %5.4f %5.4f\n",(Nodes.find(*it2)-> second)->Coord.at(0), (Nodes.find(*it2)-> second)->Coord.at(1), double(0) );
		fprintf (fOut, "%5.4f %5.4f %5.4f\n",(Nodes.at( CadhPack.List.at(*it2).firstAdhPoint ))->Coord.at(0), (Nodes.at( CadhPack.List.at(*it2).firstAdhPoint ))->Coord.at(1), double(0) );
		fprintf (fOut, "%5.4f %5.4f %5.4f\n",(Nodes.at( CadhPack.List.at(*it2).secondAdhPoint ))->Coord.at(0), (Nodes.at( CadhPack.List.at(*it2).secondAdhPoint ))->Coord.at(1), double(0) );
	}

	fprintf (fOut, "\n");
	fprintf (fOut, "LINES %d %d\n",CadhPack.ActiveCadh.size(), 3*CadhPack.ActiveCadh.size());
  	for (c = 0; c < CadhPack.ActiveCadh.size(); c++) 
    fprintf(fOut, "%d %d %d \n", 2, 2*c,2*c+1);
		
	fprintf (fOut, "\n");	
//	fprintf (fOut, "POINT_DATA %d\n",CadhPack.ActiveCadh.size()*2);
//	fprintf (fOut, "scalars Bound int\n");
//	fprintf (fOut, "LOOKUP_TABLE default\n");

	double auxX,auxY;


	fprintf (fOut, "POINT_DATA %d \n", CadhPack.ActiveCadh.size()*2);
	fprintf (fOut, "VECTORS Force_cadh double\n");
	for (j=0; j < CadhPack.ActiveCadh.size(); ++j)  
	{
		IDAux=CadhPack.ActiveCadh.at(j);
		NodesVDir(CadhPack.List.at(IDAux).firstAdhPoint ,CadhPack.List.at(IDAux).secondAdhPoint, vDir);

		auxX=CadhPack.List.at(IDAux).force*vDir.at(0);

		auxY=CadhPack.List.at(IDAux).force*vDir.at(1);
		

		fprintf (fOut, "%4.6f %4.6f 0.0\n", auxX ,auxY );

		fprintf (fOut, "%4.6f %4.6f 0.0\n", auxX*double(-1) ,auxY*double(-1) );

	}

	fprintf (fOut, "\n");

	fprintf (fOut, "CELL_DATA %d \n", CadhPack.ActiveCadh.size());
	fprintf (fOut, "SCALARS force_(pN) double\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	for (j=0; j < CadhPack.ActiveCadh.size(); ++j)  
	{
		IDAux=CadhPack.ActiveCadh.at(j);
		fprintf (fOut, "%5.4f ", CadhPack.List.at(IDAux).force*1000.);
	}
	fprintf (fOut, "\n");
	fprintf (fOut, "\n");	
//	fprintf (fOut, "CELL_DATA %d \n", CadhPack.ActiveCadh.size());
	fprintf (fOut, "SCALARS ubProb double\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	for (j=0; j < CadhPack.ActiveCadh.size(); ++j)  
	{
		IDAux=CadhPack.ActiveCadh.at(j);
		fprintf (fOut, "%5.4f ", CadhPack.List.at(IDAux).ubProb);
	}

	fprintf (fOut, "\n");
	fprintf (fOut, "\n");	
//	fprintf (fOut, "CELL_DATA %d \n", CadhPack.ActiveCadh.size());
	fprintf (fOut, "SCALARS Cluster_number int\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	for (j=0; j < CadhPack.ActiveCadh.size(); ++j)  
	{
		IDAux=CadhPack.ActiveCadh.at(j);
		fprintf (fOut, "%d ", CadhPack.List.at(IDAux).cadherinNumber);
	}
	fclose(fOut);



	
	char pref3[]="Cadherin_";	
		
	//Imprime archivo vtk
  	sprintf(name, "%s%d%s", pref3, ctParaview , extension);
 	fOut = fopen(name,"w");   //para reescribir se usa "a"
	fprintf (fOut, "# vtk DataFile Version 4.1\n");
	fprintf (fOut, "Cadherin data\n");
	fprintf (fOut, "ASCII\n");
	fprintf (fOut, "DATASET POLYDATA\n");
	fprintf (fOut, "POINTS %d double\n",CadhPack.List.size()*2);
	
	for (auto it2=0; it2 < CadhPack.List.size(); ++it2) 	
	{
		//fprintf (fOut, "%5.4f %5.4f %5.4f\n",(Nodes.find(*it2)-> second)->Coord.at(0), (Nodes.find(*it2)-> second)->Coord.at(1), double(0) );
		fprintf (fOut, "%5.4f %5.4f %5.4f\n",(Nodes.at( CadhPack.List.at(it2).firstAdhPoint ))->Coord.at(0),
											 (Nodes.at( CadhPack.List.at(it2).firstAdhPoint ))->Coord.at(1), double(0) );
		fprintf (fOut, "%5.4f %5.4f %5.4f\n",(Nodes.at( CadhPack.List.at(it2).secondAdhPoint ))->Coord.at(0), 
											 (Nodes.at( CadhPack.List.at(it2).secondAdhPoint ))->Coord.at(1), double(0) );
	}

	fprintf (fOut, "\n");
	fprintf (fOut, "LINES %d %d\n",CadhPack.List.size(), 3*CadhPack.List.size());
  	for (c = 0; c < CadhPack.List.size(); c++) 
    fprintf(fOut, "%d %d %d \n", 2, 2*c,2*c+1);
		
	fprintf (fOut, "\n");
	fprintf (fOut, "POINT_DATA %d\n",CadhPack.List.size()*2);
	fprintf (fOut, "scalars ID int\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	
	for (c = 0; c < CadhPack.List.size(); c++)  
	{
		fprintf (fOut, "%d ", 2*c);
		fprintf (fOut, "%d ", 2*c+1);
		
	}

	
	fprintf (fOut, "\n");	
	fprintf (fOut, "CELL_DATA %d \n", CadhPack.List.size());
	fprintf (fOut, "SCALARS bound double\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	for (j=0; j < CadhPack.List.size(); ++j)  fprintf (fOut, "%d ", CadhPack.List.at(j).active);
	
	fclose(fOut);

	

	// Print GAPS
	sprintf(name, "GAPS_%d%s", ctParaview , extension);
 	fOut = fopen(name,"w");   //para reescribir se usa "a"
	fprintf (fOut, "# vtk DataFile Version 4.1\n");
	fprintf (fOut, "Cell data\n");
	fprintf (fOut, "ASCII\n");
	fprintf (fOut, "DATASET POLYDATA\n");
	int i, IDGap;
	int ID;
	int ct=0;
	int sizeGap=0; 
	for (i=0;i< ActiveGAP.size();++i)
	{
		IDGap=ActiveGAP[i];
		sizeGap+=GAP.at(IDGap).CurrentShapeIDList.size();
	}

	fprintf (fOut, "POINTS %d double\n",sizeGap);
	for (i=0;i< ActiveGAP.size();++i)
	{
		IDGap=ActiveGAP[i];
		for (j=0; j<GAP.at(IDGap).CurrentShapeIDList.size(); j++) 
		{
			ID=GAP.at(IDGap).CurrentShapeIDList[j];
			fprintf (fOut, "%5.4f %5.4f %5.4f\n", 
					  Nodes[ID]->Coord.at(0), Nodes[ID]->Coord.at(1), double(0) );
		}
	}

	fprintf (fOut, "\n");
	fprintf (fOut, "POLYGONS %d %d\n",ActiveGAP.size(), ActiveGAP.size()+sizeGap);

	for (i=0;i< ActiveGAP.size();++i)
	{
		IDGap=ActiveGAP[i];
		fprintf(fOut, "%d ",GAP.at(IDGap).CurrentShapeIDList.size());
		for (j=0; j<GAP.at(IDGap).CurrentShapeIDList.size(); j++) 
		{  
    		fprintf(fOut, "%d ", ct );
    		ct++;
    	}
    	fprintf (fOut, "\n");
    }
    fprintf (fOut, "\n");

    fprintf (fOut, "POINT_DATA %d\n",sizeGap);
	fprintf (fOut, "scalars ID int\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	
	for (i=0;i< ActiveGAP.size();++i)
	{
		IDGap=ActiveGAP[i];
		for (j=0; j<GAP.at(IDGap).CurrentShapeIDList.size(); j++) 
		{
			ID=GAP.at(IDGap).CurrentShapeIDList[j];
			fprintf (fOut, "%d ",ID); 
					  
		}
	}
	
	fprintf (fOut, "\n");	

	fprintf (fOut, "CELL_DATA %d \n",ActiveGAP.size());
	fprintf (fOut, "SCALARS Area double\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	for (i=0;i<ActiveGAP.size(); ++i)
	{
		IDGap=ActiveGAP[i];	
		fprintf (fOut, "%4.4f ", GAP.at(IDGap).Area);
	}
	fprintf (fOut, "\n");

	fprintf (fOut, "SCALARS Vertex double\n");
	fprintf (fOut, "LOOKUP_TABLE default\n");
	for (i=0;i<ActiveGAP.size();++i)
	{
		IDGap=ActiveGAP[i];
		fprintf (fOut, "%d ", int(GAP.at(IDGap).vertex));
	 }
	fprintf (fOut, "\n");
	fclose(fOut);
	



	
  }	


  void Monolayer::print_data_intern(SpMat K_matrix_cc, Eigen::VectorXd &RNodes, Eigen::VectorXd &FInt , 
  								Eigen::VectorXd &Damp, Eigen::VectorXd &XIncTotal,	Eigen::VectorXd &XIncIt,int divider, int iter)
{

    FILE*fOut;
	char name[255];
	CheckNodesIncrement();
	sprintf(name,"ayuda_intern.txt");
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	fprintf (fOut," STEP: %d ******************************\n", step); 
	fprintf (fOut," ITER: %d \n", iter); 
	fprintf (fOut," DIV: %d \n", divider); 
	fprintf (fOut," \n");
	fprintf (fOut,"K_MATRIX \n");
	for (int k=0; k<K_matrix_cc.outerSize(); ++k)
  		for (SpMat::InnerIterator it(K_matrix_cc,k); it; ++it)
  		{
  		   	fprintf (fOut,"  K(i+ %d, %d +1) =   %4.8f \n", it.row(),it.col(),it.value());
		}


	fprintf (fOut," \n");
	fprintf (fOut,"Nudo     Residuo=(Stress por def)+   		Total inc de X     		-DAmper effect	    	 Incremento de X   Intern force  		 F_ext \n");
	
  	for (int i=0; i< RNodes.size() ; ++i )
		fprintf (fOut," [%4d]\t%4.8f\t%4.8f\t		%4.8f		%4.8f           %4.8f          F(i+ %4.8f )\n", 
					i, RNodes(i),XIncTotal(i), Damp(i) ,XIncIt(i),  FInt(i), ForceStep(i) );
 
	fprintf (fOut," \n");
	fprintf (fOut," R_mod: %4.8f \n", RNodes.dot(RNodes));
	fprintf (fOut," \n");

   	std::cout << std::endl;
	fclose(fOut);

}


  void Monolayer::print_data_intern( Eigen::VectorXd &FInt , Eigen::VectorXd &CInt, 
  									Eigen::VectorXd &XFext,	Eigen::VectorXd &Xvisc, int iter)
{

    FILE*fOut;
	char name[255];
	CheckNodesIncrement();
	sprintf(name,"ayuda_intern.txt");
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	fprintf (fOut," STEP: %d ******************************\n", step); 
	fprintf (fOut," ITER: %d \n", iter); 


	fprintf (fOut," \n");
	fprintf (fOut,"Nudo\tF_int\t\tC_int\t\tX_f\t\tX_visc\t\tX_stepS\t\tF_ext \n");
	
  	for (int i=0; i< FInt.size() ; ++i )
		fprintf (fOut,"[%4d]\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t%4.6f \n", 
					i, FInt(i),CInt(i), XFext(i) ,Xvisc(i),  XIncLastStep(i), ForceStep(i) );
 
 
	fclose(fOut);

}


void Monolayer::print_data_intern2(Eigen::VectorXd &FInt,Eigen::VectorXd &Faux , Eigen::VectorXd &Fvisc,
  						Eigen::VectorXd &XFext, Eigen::VectorXd &Xvisc, std::vector<double> &XIncIter, 
  						Eigen::VectorXd &XIncLastStep,int iter, int totalIter, Eigen::VectorXd &CInt)
{

    FILE*fOut;
	char name[255];
	CheckNodesIncrement();
	sprintf(name,"ayuda_intern.txt");
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	fprintf (fOut," STEP: %d \n",step);  
	fprintf (fOut," STEP : %d of %d******************************\n", iter, totalIter);
	


	fprintf (fOut," \n");
	fprintf (fOut,"Nudo\tF_int\t\tF_step\t\tF_tot\t\tX_tot\t\tCInt\t\tFvisc\t\tXvisc\t\tX_final\t\tX_acc_step\t\tXIncTotalCheck \n");
	
  	for (int i=0; i< FInt.size() ; ++i )
		fprintf (fOut,"[%4d]\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t%4.6f \n", 
					i, FInt(i),ForceStep(i), Faux(i) ,XFext(i), CInt(i), Fvisc(i), Xvisc(i),XIncIter[i], XIncLastStep(i), XIncTotalCheck(i)  );
 
 
	fclose(fOut);

}

void Monolayer::print_data_intern2(Eigen::VectorXd &FInt,Eigen::VectorXd &Faux , Eigen::VectorXd &Fvisc,
  						Eigen::VectorXd &XFext, Eigen::VectorXd &Xvisc, std::vector<double> &XIncIter, Eigen::VectorXd &CInt)
{

    FILE*fOut;
	char name[255];
	CheckNodesIncrement();
	sprintf(name,"ayuda_intern.txt");
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	fprintf (fOut," STEP: %d ******************************\n", step); 
	


	fprintf (fOut," \n");
	fprintf (fOut,"Nudo\tF_int\t\tF_ext\t\tF_tot\t\tX_tot\t\tCInt\t\tFvisc\t\tXvisc\t\tX_final\t\tXIncTotalCheck \n");
	
  	for (int i=0; i< FInt.size() ; ++i )
		fprintf (fOut,"[%4d]\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t%4.6f\t%4.6f\n", 
					i, FInt(i),ForceStep(i), Faux(i) ,XFext(i), CInt(i), Fvisc(i), Xvisc(i),XIncIter[i], XIncTotalCheck(i));
 
 
	fclose(fOut);

}


 void Monolayer::print_data_intern2( Eigen::VectorXd &FNew )
{

    FILE*fOut;
	char name[255];
	CheckNodesIncrement();
	sprintf(name,"ayuda_intern.txt");
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	fprintf (fOut," STEP: %d ******************************\n", step); 
	fprintf (fOut," GENERATION ******************************\n"); 

	fprintf (fOut," \n");
	fprintf (fOut,"Nudo\tOld_Force\tNew_force\tIncrementForce \n");
	
  	for (int i=0; i< FNew	.size() ; ++i )
		fprintf (fOut,"[%4d]\t%4.6f\t%4.6f\t%4.6f \n", 
					i, ForceLastStep(i),FNew(i), ForceInc(i)  );
 
 
	fclose(fOut);

}


void Monolayer::matlabPlot(char *index_cluster)
{
	char orden[1024];        

    sprintf(orden, "matlab -nodisplay -nosplash -r 'plotGaps(%s, %4.6f, %d,%4.6f )'", index_cluster, TimeStep, TimeSlots, GapMinTime); 

    //sprintf(orden, "python plotGaps.py"); 
    //sprintf(orden, "matlab -nodisplay -nosplash -r 'gaps_cluster(%s)'", index_cluster); 
    system (orden); 

}




void Monolayer::print_data_intern(Eigen::MatrixXd &Amat, SpMat K_matrix_cc, Eigen::VectorXd &RNodes, Eigen::VectorXd &FInt , 
  									Eigen::VectorXd &Damp, Eigen::VectorXd &XIncTotal,	Eigen::VectorXd &XIncIt, int iter)
{

    FILE*fOut;
	char name[255];
	CheckNodesIncrement();
	
	sprintf(name,"ayuda_intern.txt");
	fOut = fopen(name,"a");   //para reescribir se usa "a"
	fprintf (fOut," STEP: %d ******************************\n", step); 
	fprintf (fOut," ITER: %d \n", iter); 
	fprintf (fOut," \n");
	fprintf (fOut,"K_MATRIX \n");
	for (int k=0; k<K_matrix_cc.outerSize(); ++k)
  		for (SpMat::InnerIterator it(K_matrix_cc,k); it; ++it)
  		{
  		   	fprintf (fOut,"  K(i+ %d, %d +1) =   %4.8f \n", it.row(),it.col(),it.value());
		}

	fprintf (fOut," \n");
	fprintf (fOut,"A_MATRIX \n");
	for (int k=0; k<(2* meshSize); ++k)
	{
  		for (int i=0; i<(2* meshSize); ++i)
  		{
  		   	fprintf (fOut,"  K(1+ %d, %d +1) =   %4.8f \n", k,i,Amat(k,i));
		}
	}

	Eigen::MatrixXd InvA;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(Amat);
//   if (eigensolver.info() != Eigen::Success) abort();
   std::cout << "The eigenvalues of A are: " << eigensolver.eigenvalues() << std::endl;
   std::cout << "The determinant of A is " << Amat.determinant() << std::endl;
   fprintf (fOut,"\n  The determinant is:  %4.8f \n", Amat.determinant());

   	InvA= Amat.inverse();

	
   	Damp=InvA*RNodes;

   	fprintf (fOut," \n");
	fprintf (fOut,"Nudo       incX_con inverse		 F_ext \n");
	
  	for (int i=0; i< RNodes.size() ; ++i )
		fprintf (fOut,"  %d        %4.4f    		%4.4f\n", 
					i,  Damp(i), ForceStep(i) );
 
	fprintf (fOut," \n");
	fprintf (fOut," R_mod: %4.8f \n", RNodes.dot(RNodes));
	fprintf (fOut," \n");

   


	fclose(fOut);

}


/*



////////////////////////////////////////

// Apply boundary conditions

	 cellCt=0;
	 row=nCellRow;
	 int encastredID, encastredIDmax;
	 // Center and upper zone
	 for (j=0; j < column-1; ++j)
	 {
	 	row=nCellRow-j;
	 	
	 	// bdry in x direction
		encastredID=GETVertexID(cellCt,4);
		encastredIDmax=GETVertexID(cellCt,5);
		
		for (i=encastredID; i<= encastredIDmax;++i )
		{
			BdryIDs.push_back(i*2);

		}


		cellCt=cellCt+row;

		encastredID=GETVertexID(cellCt-1,1);
		encastredIDmax=GETVertexID(cellCt-1,2);

		for (i=encastredID; i<= encastredIDmax;++i )
		{
			BdryIDs.push_back(i*2);

		}

	 }

	 // Upper cell
	 encastredID=GETVertexID(cellCt,4);
	 encastredIDmax=GETVertexID(cellCt,5);
		
	 for (i=encastredID; i<= encastredIDmax;++i )
		BdryIDs.push_back(i*2);

	 
	 encastredID=GETVertexID(cellCt,1);
	 encastredIDmax=GETVertexID(cellCt,2);

	 for (i=encastredID; i<= encastredIDmax;++i )
		BdryIDs.push_back(i*2);	

	// bdry in y direction
	 encastredID=GETVertexID(cellCt,0);
	 BdryIDs.push_back(encastredID*2+1);
	 cellCt++;

	 // Lower zone
	 row=nCellRow;
	 for (j=1; j < column-1; ++j)   // Check if it -1
	 {
	 	row=nCellRow-j;
		// bdry in x direction
		encastredID=GETVertexID(cellCt,4);
		encastredIDmax=GETVertexID(cellCt,5);
		
		for (i=encastredID; i<= encastredIDmax;++i )
		{
			BdryIDs.push_back(i*2);

		}
		cellCt=cellCt+row;

		encastredID=GETVertexID(cellCt-1,1);
		encastredIDmax=GETVertexID(cellCt-1,2);

		for (i=encastredID; i<= encastredIDmax;++i )
		{
			BdryIDs.push_back(i*2);

		}	
	 }

	 // Lower cell
	 encastredID=GETVertexID(cellCt,4);
	 encastredIDmax=GETVertexID(cellCt,5);
		
	 for (i=encastredID; i<= encastredIDmax;++i )
		BdryIDs.push_back(i*2);

	 
	 encastredID=GETVertexID(cellCt,1);
	 encastredIDmax=GETVertexID(cellCt,2);

	 for (i=encastredID; i<= encastredIDmax;++i )
		BdryIDs.push_back(i*2);	


	// bdry in y direction
	 encastredID=GETVertexID(cellCt,3);
	 BdryIDs.push_back(encastredID*2+1);
	 cellCt++;



 */





/*
void Monolayer::PrintGAP(char *index_cluster)
{
	int ctVertex=0;
	int ctInt=0;
	char name[100];
	char orden[1024]; 
	int ct=0;
	int IDGap;
	double vertexTime=double(0);
	double interfaceTime=double(0);
	FILE*fOut;

	for (auto i=0; i< InactiveGAP.size(); ++i)
	{
		IDGap=InactiveGAP[i];
		if ( GAP.at(IDGap).vertex )
		{
			sprintf(name, "GapVertex_%s_%d.txt" , index_cluster, ctVertex);			// Number of Gaps generated
			fOut = fopen(name,"a");   //para reescribir se usa "a"
			ctVertex++;
			vertexTime+=GAP.at(IDGap).Area.size()*TimeStep;

		}
		else
		{
			sprintf(name, "GapInterface_%s_%d.txt" , index_cluster, ctInt);			// Number of Gaps generated
			fOut = fopen(name,"a");   //para reescribir se usa "a"
			ctInt++;
			interfaceTime+=GAP.at(IDGap).Area.size()*TimeStep;
		}
		for (ct=0; ct < GAP.at(IDGap).Area.size(); ++ct )
			fprintf (fOut, "%d %d %d %d %4.4f %4.4f %4.4f %4.4f %4.4f\n ",GAP.at(IDGap).startStep+ct,  GAP.at(IDGap).nAdhPoint[ct], GAP.at(IDGap).InitPoint[ct], GAP.at(IDGap).LastPoint[ct], 
						GAP.at(IDGap).Area[ct], GAP.at(IDGap).Length[ct], GAP.at(IDGap).HeightAver[ct], GAP.at(IDGap).HeightMax[ct], GAP.at(IDGap).DistVertex[ct] );   // Print mean traction and last traction		    

		fclose(fOut);
		
		
	}

	for (auto i=0; i< ActiveGAP.size(); ++i)
	{
		IDGap=ActiveGAP[i];
		if ( GAP.at(IDGap).vertex )
		{
			sprintf(name, "GapVertex_%s_%d.txt" , index_cluster, ctVertex);			// Number of Gaps generated
			fOut = fopen(name,"a");   //para reescribir se usa "a"
			ctVertex++;
			vertexTime+=GAP.at(IDGap).Area.size()*TimeStep;
		}
		else
		{
			sprintf(name, "GapInterface_%s_%d.txt" , index_cluster, ctInt);			// Number of Gaps generated
			fOut = fopen(name,"a");   //para reescribir se usa "a"
			ctInt++;
			interfaceTime+=GAP.at(IDGap).Area.size()*TimeStep;			
		}
		for (ct=0; ct < GAP.at(IDGap).Area.size(); ++ct )
			fprintf (fOut, "%d %d %d %4.4f %4.4f %4.4f %4.4f %4.4f\n ",GAP.at(IDGap).startStep+ct, GAP.at(IDGap).InitPoint[ct], GAP.at(IDGap).nAdhPoint[ct], 
						GAP.at(IDGap).Area[ct], GAP.at(IDGap).Length[ct], GAP.at(IDGap).HeightAver[ct], GAP.at(IDGap).HeightMax[ct], GAP.at(IDGap).DistVertex[ct] );   // Print mean traction and last traction		    

		fclose(fOut);
		
		
	}

	std::cout << "Vertex -> Number: " << ctVertex << ";  Seconds:"<< vertexTime <<std::endl;
	std::cout << "Interface -> Number: " << ctInt << ";  Seconds:"<< interfaceTime <<std::endl;
}

*/