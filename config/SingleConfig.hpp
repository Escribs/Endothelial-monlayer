#ifndef _H_CONFIGURATION_SINGLETON_
#define _H_CONFIGURATION_SINGLETON_


// STL C++ Library
#include <iostream>
#include <fstream>
#include <vector>
#include <string>


class SingleConfig
{
    public:
        static SingleConfig& getConfig();
       double cadhLengthBal, radius, sepAdhPoint ,sepInterface ,aristLength, 
              TimeStep, cadhDensity, cadhEA, cadhVisc, KonRate, bindLimitDistance, mediumViscCoef;
          
       double bdryEA, membEA, membVisc, bdryVisc;

       double maxCtForce, maxRadForce, maxFlowForce, dirFlowForce, stepForce, stepSmoothForce, maximunDispAllowed;

       double KCatch, phiS, phiC, fAp;  

       double KBendingCadh, KBendingMemb;

       double gapThreshold; 

       int nArists,nCell, TimeSlots, nAdhPerPoint;

       bool cellCenter, testFormulation, printFilesVTK;

       int smallDefResolution, plotFreq;

       int clusterMaxNumber, initialCadherinNumber;

       int extraBdry;

       int nPropagationMemb, nPropagationBdry; 
       double ratioZeroRad, ratioZeroCt;  


       // Remodeling

       double cellExtraGrowth, Kremodel, remodelRatioVel, remodelRatioVelVS;
       bool remodelUB;

       // Protrusion

        double maxProtForce, ratioZeroProt, stepForceProtr;
        int nPropagationProt;



       
       bool readFromFile();


    private:
        SingleConfig();                   // Constructor? (the {} brackets) are needed here.

       // Configuration file to read
        static std::string inputFile;

        // C++ 11
        // =======
        // We can use the better technique of deleting the methods
        // we don't want.
        //SingleConfig(SingleConfig const&)               = delete;
       // void operator=(SingleConfig const&)             = delete;

        // variables 

};
#endif 