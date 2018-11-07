// ##################################################
// #   header.hpp - finally revised on Jun 2018     #
// #   coded by Jorge Escribano                     #
// #   Copyright (C) 2016-2018, Jorge Escribano,    #
// #   All rights reserved.                         #
// ##################################################


#ifndef HEADER_HPP
#define HEADER_HPP

 #include <ctime>  //For the seed of the aleatory builder
 #include <cmath> // For trigonometric functions
 #include <iomanip>   // for using  <<  and rand() 
 //#include <time.h>  //For the seed of the aletory builder
 #include <iostream> //For std
 #include <math.h>
 #include <stdlib.h> //For rand()	
 #include <fstream> //For fclose
 #include <vector>
 #include <stdio.h>	//strcpy
 #include <string>	//strcpy <string.h>
 #include <algorithm>    // std::max
 //#include <boost/lexical_cast.hpp>   //Libreria locura
// Eigen includes
#include <chrono>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <memory>
#include <map>
#include <Eigen/IterativeLinearSolvers>
// For CGal
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <set>

 const double PI=3.1415926535897932;

 typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
 typedef Eigen::Triplet<double> T;

// For CGal
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 PointCGal;
typedef CGAL::Polygon_2<K> Polygon_2;

#endif /* HEADER_HPP */