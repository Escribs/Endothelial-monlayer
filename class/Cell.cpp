// ##################################################
// #   Cell.cpp - finally revised on Jun 2018     	#
// #   coded by Jorge Escribano                     #
// #   Copyright (C) 2016-2018, Jorge Escribano,    #
// #   All rights reserved.                         #
// ##################################################

// Contains the cell class

#include "Cell.hpp"


Cell::Cell(std::vector<std::vector<int>>& ConnectivityMemb, std::vector<std::vector<int>>& ConnectivityBdry, int IDCtr)
 {
 	int IDNodes, IDNodesLast, i, IDNodesAux, vertexPlusCadh;
	
	 ID=IDCtr;
	// Read variables
	 SingleConfig DataConfig=SingleConfig::getConfig();
	 cellCenter=DataConfig.cellCenter;
	 nArists=DataConfig.nArists;
	 sepAdhPoint=DataConfig.sepAdhPoint;
	 cadhLength=DataConfig.cadhLengthBal;
	 sepInterface=DataConfig.sepInterface;
	 radius=DataConfig.radius;
	 aristLength=DataConfig.aristLength;
	 int nAdhPerPoint=DataConfig.nAdhPerPoint;
	 if (nAdhPerPoint==1)
	 	vertexPlusCadh=1; // To add in vertex two cadh to each cell
	 else
	 	vertexPlusCadh=0; // To add an ad
	 

	 int nCell= DataConfig.nCell;

	 if(nCell==2)
 		nArists=1;
 	 else if(nCell==3)
 		nArists=2;


	  if (DataConfig.testFormulation)
	  {

	 	nCell=1;
	 	sepAdhPoint=1;
	 	nArists=1;
	 	aristLength=1;

	 	std::vector<double> Coord(2);

	 	
	 	Coord.at(0)=double(0);
	 	Coord.at(1)=double(0);
	 	IDNodes=0;
	 	PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint)) );
	 	Coord.at(0)=double(0);
	 	Coord.at(1)=double(1);
	 	IDNodes=1;
	 	PointBdry.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint)) );
	 	Coord.at(0)=double(1);
	 	Coord.at(1)=double(0);
	 	IDNodes=2;
	 	PointBdry.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint)) );
	 

	 	ConnectivityBdry[0].push_back(2);	
		ConnectivityBdry[1].push_back(0);
		ConnectivityBdry[0].push_back(1);	
		ConnectivityBdry[1].push_back(0);
	

/*	
	 	Coord.at(0)=double(1);
	 	Coord.at(1)=double(0);
	 	IDNodes=1;
	 	PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes)) );

	 	Coord.at(0)=double(2);
	 	Coord.at(1)=double(-1);
	 	IDNodes=4;
	 	PointBdry.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes)) );
	 	Coord.at(0)=double(2);
	 	Coord.at(1)=double(1);
	 	IDNodes=5;
	 	PointBdry.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes)) );

	 	ConnectivityBdry[0].push_back(4);	
		ConnectivityBdry[1].push_back(1);
		ConnectivityBdry[0].push_back(5);	
		ConnectivityBdry[1].push_back(1);

		ConnectivityMemb[0].push_back(0);
		ConnectivityMemb[1].push_back(1);	
//****************

		Coord.at(0)=double(0);
	 	Coord.at(1)=double(0);
	 	IDNodes=0;
	 	PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes)) );
	 	Coord.at(0)=double(1);
	 	Coord.at(1)=double(0);
	 	IDNodes=1;
	 	PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes)) );
	 	Coord.at(0)=double(0);
	 	Coord.at(1)=double(2);
	 	IDNodes=2;
	 	PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes)) );
	 	Coord.at(0)=double(1);
	 	Coord.at(1)=double(2);
	 	IDNodes=3;
	 	PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes)) );

	 	Coord.at(0)=double(-1);
	 	Coord.at(1)=double(1);
	 	IDNodes=4;
	 	PointBdry.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes)) );
	 	Coord.at(0)=double(2);
	 	Coord.at(1)=double(1);
	 	IDNodes=5;
	 	PointBdry.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes)) );

	 	ConnectivityBdry[0].push_back(4);	
		ConnectivityBdry[1].push_back(0);
		ConnectivityBdry[0].push_back(4);	
		ConnectivityBdry[1].push_back(2);
		ConnectivityBdry[0].push_back(5);	
		ConnectivityBdry[1].push_back(1);
		ConnectivityBdry[0].push_back(5);	
		ConnectivityBdry[1].push_back(3);

		ConnectivityMemb[0].push_back(0);
		ConnectivityMemb[1].push_back(1);	
		ConnectivityMemb[0].push_back(2);
		ConnectivityMemb[1].push_back(3);	
*/
	 }

	 else

	 {
	 	std::vector<int> MembPointsID;

		// Number of adhesion points at each arist  
		 nAdhPoint=floor(aristLength/sepAdhPoint)+1; 
		 aristLength=double(nAdhPoint*sepAdhPoint);
	
	    // Create adhesion points
		 std::vector<double> Coord(2);
		 std::vector<double> CoordAux(2);
		 std::vector<double> VDir(2);
		
	
		 PointsDirection(Coord,VDir,true);
		 
	
	//*********************************// Tener en cuenta que cuando se pongan dos aristas hay un punto que se repite habra que poner un  -1 en las cuentas
		// Create Membrane points and their connectivity 
	 
		 i=0;

		 if (ID==3)
		 	IDNodes=nAdhPoint*2;
		 else
			 IDNodes=2*i+(ID-1);				// How I have designed the mesh. Look notebook

	 	 PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint)) );
	 	 MembPointsID.push_back(IDNodes);
		 IDNodesLast=IDNodes;	 
	
		 for (i=1; i<nAdhPoint; ++i)
		 {
		 	if (ID==3)
		 		IDNodes=(nAdhPoint*2-1)+4*i;	// Goes with cell 1 interface
		 	else
			 	IDNodes=2*i+(ID-1);				// How I have designed the mesh. Look notebook

		 	Coord.at(0)+=VDir.at(0)*sepAdhPoint;
			Coord.at(1)+=VDir.at(1)*sepAdhPoint;
			
			ConnectivityMemb[0].push_back(IDNodesLast);	
			ConnectivityMemb[1].push_back(IDNodes);
			PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint)) );
			MembPointsID.push_back(IDNodes);
			IDNodesLast=IDNodes;
		 } 
	
		 // Create a different one for a third cell
		 if (nArists==2)
		 {
		 	// Fill after first tests
		 	PointsDirection3cells(VDir);
		 	int RefNode;
		 	
		 	if (ID==3)
		 	{
		 		RefNode=nAdhPoint*2;
		 		Coord.at(0)=sepInterface/2.;    
				Coord.at(1)=sepAdhPoint*(nAdhPoint-1)+sepInterface*sin(3.1416/3);  
				IDNodesLast=  RefNode;
		 	}
		 	else 
		 		RefNode=IDNodesLast-1;
		 	
		 	for (i=1; i<nAdhPoint; ++i)
		 	{
		 		IDNodes=RefNode+4*i;				// How I have designed the mesh. Look notebook
		 		Coord.at(0)+=VDir.at(0)*sepAdhPoint;
				Coord.at(1)+=VDir.at(1)*sepAdhPoint;
			
				ConnectivityMemb[0].push_back(IDNodesLast);	
				ConnectivityMemb[1].push_back(IDNodes);
				PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint)) );
				MembPointsID.push_back(IDNodes);
				IDNodesLast=IDNodes;
			}
			CoordAux=Coord; 
			IDNodesAux=IDNodesLast;
		 }
	 	
		 //
		 
	 	 PointsDirection(Coord,VDir,false); // To reset Coord
		 if (cellCenter) // Create nodes in radial configuration
		 {
		 	if (nCell==2)
		 	{
		 		IDNodesLast=(ID-1)+nArists*nCell*nAdhPoint;
		 		Coord.at(1)+=VDir.at(1)*sepAdhPoint*(nAdhPoint-1)/double(2);
		 	}
		 	else  //nCell=3
		 	{
		 		IDNodesLast=(ID-1)+nCell*(nAdhPoint*2-1);

		 		if (ID==3)
		 			Coord.at(1)+=radius;
				else if (ID==2)
				{
					Coord.at(0)=radius*cos(3.1416/6.)+sepInterface;
					Coord.at(1)=sepAdhPoint*(nAdhPoint-1)-radius*sin(3.1416/6.);
				}
				else if (ID==1)
				{
					Coord.at(0)=-radius*cos(3.1416/6.);
					Coord.at(1)=sepAdhPoint*(nAdhPoint-1)-radius*sin(3.1416/6.);
				}
			}

			PointBdry.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodesLast, nAdhPerPoint)) );
			

			for (i=0; i<MembPointsID.size(); ++i)
		 	{	
		 		//IDNodes=2*i+(ID-1);
		 		IDNodes=MembPointsID.at(i);
				ConnectivityBdry[0].push_back(IDNodesLast);	
				ConnectivityBdry[1].push_back(IDNodes);
			}
		 }
	
		 else  // Create nodes in wall configuration
		 {
	 		i=0;
		 	IDNodes=2*i+(ID-1)+nArists*nCell*nAdhPoint;				// How I have designed the mesh. Look notebook
	 	 	PointBdry.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint)) );
		 	ConnectivityBdry[0].push_back(IDNodes);	
			ConnectivityBdry[1].push_back(IDNodes-nArists*nCell*nAdhPoint);	 
	
		 	for (i=1; i<nAdhPoint; ++i)
		 	{	
	
				IDNodes=2*i+(ID-1)+nArists*nCell*nAdhPoint;				// How I have designed the mesh. Look notebook
		 		Coord.at(0)+=VDir.at(0)*sepAdhPoint;
				Coord.at(1)+=VDir.at(1)*sepAdhPoint;
				
				ConnectivityBdry[0].push_back(IDNodes);	
				ConnectivityBdry[1].push_back(IDNodes-nArists*nCell*nAdhPoint);
				PointBdry.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint)) );
		 	}
		}

		// Add boundary in y direction so K matrix is not singular 
		
		PointsDirection(Coord,VDir,true); // To reset Coord
		
	 	if (ID==3)
		{
			i=nAdhPoint-1;
			IDNodes=(nAdhPoint*2-1)+4*i;
	 		Coord.at(0)+=VDir.at(0)*sepAdhPoint*nAdhPoint;
			Coord.at(1)+=VDir.at(1)*sepAdhPoint*nAdhPoint;
		}
		else
		{
			i=0;
			IDNodes=2*i+(ID-1);
			Coord.at(0)+=VDir.at(0);
			Coord.at(1)-=VDir.at(1)*sepAdhPoint;
		}

		ConnectivityBdry[1].push_back(IDNodes);	  // Memb point

		if (cellCenter)
		{
			IDNodes=(ID-1)+nArists*nCell*nAdhPoint+2;
			if (nCell==3)
				IDNodes=(nCell*nAdhPoint*2)+(ID-1);

		}
		else
		{
			IDNodes=(ID-1)+nArists*nCell*nAdhPoint*2;
		}

		ConnectivityBdry[0].push_back(IDNodes);   // Bdry point
		PointBdry.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint)) );
		/*
		PointsDirection(Coord,VDir,true); // To reset Coord
		i=nAdhPoint-1;  
		IDNodes=2*i+(ID-1);	// Last up node
		ConnectivityBdry[1].push_back(IDNodes);	  // Memb point
		Coord.at(0)+=VDir.at(0);
		Coord.at(1)+=VDir.at(1)*sepAdhPoint*(nAdhPoint);
		if (cellCenter)
		{
			IDNodes=(ID-1)+nArists*nCell*nAdhPoint+4; 
		}
		else
		{
			i=nAdhPoint+1;
			IDNodes=2*i+(ID-1)+nArists*nCell*nAdhPoint;
		}
		ConnectivityBdry[0].push_back(IDNodes);   // Bdry point
		PointBdry.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes)) ); 
		*/
	

		if (nCell==3)
		{
		 	PointsDirection3cells(VDir);			
			Coord.at(0)=CoordAux.at(0)+VDir.at(0)*sepAdhPoint;
			Coord.at(1)=CoordAux.at(1)+VDir.at(1)*sepAdhPoint;
			ConnectivityBdry[1].push_back(IDNodesAux);	  // Memb point
			IDNodes=(nCell*nAdhPoint*2)+(ID-1)+3;
			ConnectivityBdry[0].push_back(IDNodes);	  // Bdry point
			PointBdry.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint)) );		
		}

	}

	POINTSCHECK();
 }

 void Cell::PointsDirection(std::vector<double> &Coord, std::vector<double> &VDir, bool Membrane)	
 {
 	if (Membrane)
 	 {	
 	 	if (ID==1)
		{	 
			Coord.at(0)=double(0);
			Coord.at(1)=double(0);
			VDir.at(0)=double(0);
			VDir.at(1)=double(+1);
		} 
	
		else if (ID==2)
		 {
		 	Coord.at(0)=sepInterface;    
			Coord.at(1)=double(0);
			VDir.at(0)=double(0);
			VDir.at(1)=double(+1);
		 }
	
		 else if (ID==3)
		 {

		 	Coord.at(0)=sepInterface/2.;    
			Coord.at(1)=sepAdhPoint*(nAdhPoint-1)+sepInterface*sin(3.1416/3);  
		 	VDir.at(0)=-cos(3.1416/6);
			VDir.at(1)=sin(3.1416/6);
		 }
	 }
	
	else // Boundary point
	 {
		if (ID==1)
		{	 
			Coord.at(0)=double(-radius);
			Coord.at(1)=double(0);
			VDir.at(0)=double(0);
			VDir.at(1)=double(1);
		} 
	
		else if (ID==2)
		 {
		 	Coord.at(0)=sepInterface+radius;
			Coord.at(1)=double(0);
			VDir.at(0)=double(0);
			VDir.at(1)=double(+1);
		 }
	
		 else if (ID==3)
		 {
		 	Coord.at(0)=sepInterface/2.;    
			Coord.at(1)=sepAdhPoint*(nAdhPoint-1)+sepInterface*sin(3.1416/3);  
		 }
	}



 }

 void Cell::PointsDirection3cells( std::vector<double> &VDir)	
 {
 	if (ID==1)
	{	 
		VDir.at(0)=-cos(3.1416/6);
		VDir.at(1)=sin(3.1416/6);
	} 
	
	else if (ID==2)
	{
		VDir.at(0)=cos(3.1416/6);
		VDir.at(1)=sin(3.1416/6);
	}
	
	else if (ID==3)
	{

		VDir.at(0)=cos(3.1416/6);
		VDir.at(1)=sin(3.1416/6);
	}
 }




void Cell::POINTSCHECK()
  {
	int j,c, nConnectivity;
	int IDAux;
	char name[50];
	char pref[]="Cell_check";	
	char extension[]=".txt";

	 std::vector<std::shared_ptr<MemAdhPoint>>::iterator it;

	

	FILE*fOut;
	//Imprime archivo vtk
  	sprintf(name, "%s%s", pref , extension);
 	fOut = fopen(name,"a");   //para reescribir se usa "a"
	
	fprintf (fOut, "Cell: %d \n", ID);
	
	fprintf (fOut, "POINTS MEMB double\n");
	
	for (it=PointMemb.begin(); it !=PointMemb.end(); ++it) 	fprintf (fOut, "ID:%d   Coord: %5.4f %5.4f %5.4f\n", 
													(*it)->get_IDNode(), (*it)->Coord.at(0), (*it)->Coord.at(1), double(0) );

	fprintf (fOut, "POINTS Bdry double\n");
	
	for (it=PointBdry.begin(); it !=PointBdry.end(); ++it) 	fprintf (fOut, "ID:%d   Coord: %5.4f %5.4f %5.4f\n", 
													(*it)->get_IDNode(), (*it)->Coord.at(0), (*it)->Coord.at(1), double(0) );


	
	fclose(fOut);

}


Cell::Cell( std::vector<double> Center ,std::vector<std::vector<int>>& ConnectivityMemb, std::vector<std::vector<int>>& ConnectivityBdry, int IDCtr)
 {
 	int IDNodes, IDNodesLast, i, IDNodesAux, nAdhPerPoint, vertexPlusCadh;
 	std::vector<double> Coord(2);
	
	 ID=IDCtr; // Cell ID
	// Read variables
	SingleConfig DataConfig=SingleConfig::getConfig();
	cellCenter=DataConfig.cellCenter;
	nArists=DataConfig.nArists;
	sepAdhPoint=DataConfig.sepAdhPoint;
	cadhLength=DataConfig.cadhLengthBal;
	sepInterface=DataConfig.sepInterface;
	radius=DataConfig.radius;
	aristLength=DataConfig.aristLength;
	nAdhPerPoint=DataConfig.nAdhPerPoint;  // Number of cadh permited per AdhPoint
	int nCell= DataConfig.nCell;

	 if (nAdhPerPoint==1)
	 	vertexPlusCadh=1; // To add in vertex two cadh to each cell
	 else
	 	vertexPlusCadh=0; // To add an ad
	 

	nAdhPointSide=floor((aristLength+0.1)/sepAdhPoint); 
	aristLength=double(nAdhPointSide*sepAdhPoint);
	nAdhPoint=nAdhPointSide*6;

	// Cell * points in a side * 6 sides + center point 
	IDNodes=(ID)*(nAdhPoint +1);
	IDCenter=IDNodes;

	CenterPoint=Center;
	Coord=Center;

	//UnboundAdhPoint.reserve(nAdhPoint);

	//  center point
	PointBdry.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(CenterPoint, IDNodes, nAdhPerPoint)) );
	IDNodes++;
	
	// Create Faces. Works like a clock

	int auxVertexPlus;
	auxVertexPlus=vertexPlusCadh;
	Coord.at(1)+=aristLength;
	PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint+auxVertexPlus)) );
	IDNodes++;
	ConnectivityMemb[0].push_back(IDNodes-1);	
	ConnectivityMemb[1].push_back(IDNodes);

	// Cara 1
	for ( i=0; i < nAdhPointSide; i++ )
	{
		if (i==nAdhPointSide-1)
			auxVertexPlus=vertexPlusCadh;
		else
			auxVertexPlus=0;

		Coord[0]+=sepAdhPoint*cos(PI/6);
		Coord[1]-=sepAdhPoint*sin(PI/6);
		PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint+auxVertexPlus)) );
		IDNodes++;
		ConnectivityMemb[0].push_back(IDNodes-1);	
		ConnectivityMemb[1].push_back(IDNodes);
	}

	// Cara 2 
	for ( i=0; i < nAdhPointSide; i++ )
	{
		if (i==nAdhPointSide-1)
			auxVertexPlus=vertexPlusCadh;
		else
			auxVertexPlus=0;

		Coord[1]-=sepAdhPoint;
		PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint+auxVertexPlus)) );
		IDNodes++;
		ConnectivityMemb[0].push_back(IDNodes-1);	
		ConnectivityMemb[1].push_back(IDNodes);
	}

	// Cara 3
	for ( i=0; i < nAdhPointSide; i++ )
	{
		if (i==nAdhPointSide-1)
			auxVertexPlus=vertexPlusCadh;
		else
			auxVertexPlus=0;

		Coord[0]-=sepAdhPoint*cos(PI/6);
		Coord[1]-=sepAdhPoint*sin(PI/6);
		PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint+auxVertexPlus)) );
		IDNodes++;
		ConnectivityMemb[0].push_back(IDNodes-1);	
		ConnectivityMemb[1].push_back(IDNodes);
	}
	

	// Cara 4
	for ( i=0; i < nAdhPointSide; i++ )
	{
		if (i==nAdhPointSide-1)
			auxVertexPlus=vertexPlusCadh;
		else
			auxVertexPlus=0;

		Coord[0]-=sepAdhPoint*cos(PI/6);
		Coord[1]+=sepAdhPoint*sin(PI/6);
		PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint+auxVertexPlus)) );
		IDNodes++;
		ConnectivityMemb[0].push_back(IDNodes-1);	
		ConnectivityMemb[1].push_back(IDNodes);
	}
	

	// Cara 5 
	for ( i=0; i < nAdhPointSide; i++ )
	{
		if (i==nAdhPointSide-1)
			auxVertexPlus=vertexPlusCadh;
		else
			auxVertexPlus=0;

		Coord[1]+=sepAdhPoint;
		PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint+auxVertexPlus)) );
		IDNodes++;
		ConnectivityMemb[0].push_back(IDNodes-1);	
		ConnectivityMemb[1].push_back(IDNodes);
	}
	

	// Cara 6
	for ( i=0; i < nAdhPointSide-1; i++ )  
	// -1 because the first point is already added at the beggening
	{
		Coord[0]+=sepAdhPoint*cos(PI/6);
		Coord[1]+=sepAdhPoint*sin(PI/6);
		PointMemb.push_back( std::make_shared<MemAdhPoint>(MemAdhPoint(Coord, IDNodes, nAdhPerPoint)) );
		IDNodes++;
		if ( i==(nAdhPointSide-2) )  // Last one binds to the first
		{
			ConnectivityMemb[0].push_back(IDNodes-1);	
			ConnectivityMemb[1].push_back(IDCenter+1);
		}	
		else
		{
			ConnectivityMemb[0].push_back(IDNodes-1);	
			ConnectivityMemb[1].push_back(IDNodes);
		}
	}

	// Create Boundary conditions
	for (i=0; i<nAdhPoint;i++)
	{
		ConnectivityBdry[0].push_back(IDCenter);
		ConnectivityBdry[1].push_back(IDCenter+i+1);
	}
}