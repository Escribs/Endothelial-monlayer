//====================== CONFIGURATION CLASS DEFINITIONS =====================//

// HEADERS
#include "SingleConfig.hpp"
#include "PropertiesParser.hpp"



//------------------------------ STATIC VARIABLES ----------------------------//

std::string SingleConfig::inputFile = "./input.conf" ;

//------------------------------ STATIC FUNCTIONS ----------------------------//

SingleConfig::SingleConfig() 
{
    SingleConfig::readFromFile();
}

SingleConfig& SingleConfig::getConfig()
{
    static SingleConfig    instance; // Guaranteed to be destroyed.
                                  // Instantiated on first use.
    return instance;
}

//----------------------------------------------------------------------------//
//---------------------------- CONFIGURATION READER --------------------------//
//----------------------------------------------------------------------------//

bool SingleConfig::readFromFile()
{
    
     // Open input file
   	std::ifstream file( inputFile.c_str() );
    if ( ! file.is_open() ) 
    {
        std::cout   << "Error: Unable to open configuration file: " 
                    << inputFile << std::endl;
        std::abort(); 
    }


      
    //Feed properties parser with the variables to be read
    base::io::PropertiesParser prop;


    prop.registerPropertiesVar( "testFormulation",          testFormulation );

    prop.registerPropertiesVar( "printFilesVTK",            printFilesVTK );
    prop.registerPropertiesVar( "TimeSlots",                TimeSlots );
    prop.registerPropertiesVar( "TimeStep",                 TimeStep );
    prop.registerPropertiesVar( "stepForce",                stepForce);
    prop.registerPropertiesVar( "stepSmoothForce",          stepSmoothForce);   
    prop.registerPropertiesVar( "radius",                   radius);
    prop.registerPropertiesVar( "sepAdhPoint",              sepAdhPoint);
    prop.registerPropertiesVar( "sepInterface",             sepInterface);
    prop.registerPropertiesVar( "smallDefResolution",       smallDefResolution);    
    prop.registerPropertiesVar( "mediumViscCoef",           mediumViscCoef);
    prop.registerPropertiesVar( "nAdhPerPoint",             nAdhPerPoint);
    prop.registerPropertiesVar( "plotFreq",                 plotFreq);
    

    prop.registerPropertiesVar( "maximunDispAllowed",       maximunDispAllowed);
   
    prop.registerPropertiesVar( "gapThreshold",             gapThreshold);  

    prop.registerPropertiesVar( "nPropagationMemb",         nPropagationMemb);  
    prop.registerPropertiesVar( "nPropagationBdry",         nPropagationBdry);  
    prop.registerPropertiesVar( "ratioZeroRad",             ratioZeroRad);  
    prop.registerPropertiesVar( "ratioZeroCt",              ratioZeroCt);  


    prop.registerPropertiesVar( "aristLength",              aristLength );
    prop.registerPropertiesVar( "nArists",                  nArists);
    prop.registerPropertiesVar( "cellCenter",               cellCenter);
    prop.registerPropertiesVar( "nCell",                    nCell);

    prop.registerPropertiesVar( "cadhLengthBal",            cadhLengthBal );
    prop.registerPropertiesVar( "KBendingCadh",             KBendingCadh); 
    prop.registerPropertiesVar( "cadhDensity",              cadhDensity);
    prop.registerPropertiesVar( "cadhEA",                   cadhEA);
    prop.registerPropertiesVar( "cadhVisc",                 cadhVisc);
    prop.registerPropertiesVar( "KonRate",                  KonRate);
    prop.registerPropertiesVar( "clusterMaxNumber",         clusterMaxNumber);
    prop.registerPropertiesVar( "initialCadherinNumber",    initialCadherinNumber);

    prop.registerPropertiesVar( "maxCtForce",               maxCtForce);
    prop.registerPropertiesVar( "maxRadForce",              maxRadForce);

    prop.registerPropertiesVar( "maxFlowForce",             maxFlowForce);
    prop.registerPropertiesVar( "dirFlowForce",             dirFlowForce);

    prop.registerPropertiesVar( "extraBdry",                extraBdry);

    prop.registerPropertiesVar( "bindLimitDistance",        bindLimitDistance);

    prop.registerPropertiesVar( "bdryEA",                   bdryEA);
    prop.registerPropertiesVar( "membEA",                   membEA);
    prop.registerPropertiesVar( "bdryVisc",                 bdryVisc);
    prop.registerPropertiesVar( "membVisc",                 membVisc);    
    prop.registerPropertiesVar( "KBendingMemb",             KBendingMemb); 


    prop.registerPropertiesVar( "KCatch",                   KCatch);
    prop.registerPropertiesVar( "phiS",                     phiS);
    prop.registerPropertiesVar( "phiC",                     phiC);
    prop.registerPropertiesVar( "fAp",                      fAp);  




    prop.registerPropertiesVar( "cellExtraGrowth",              cellExtraGrowth);
    prop.registerPropertiesVar( "Kremodel",                     Kremodel);
    prop.registerPropertiesVar( "remodelRatioVel",              remodelRatioVel);
    prop.registerPropertiesVar( "remodelRatioVelVS",            remodelRatioVelVS);  
    prop.registerPropertiesVar( "remodelUB",                    remodelUB);  


    prop.registerPropertiesVar( "maxProtForce",              maxProtForce);
    prop.registerPropertiesVar( "ratioZeroProt",             ratioZeroProt);
    prop.registerPropertiesVar( "stepForceProtr",            stepForceProtr);
    prop.registerPropertiesVar( "nPropagationProt",          nPropagationProt);  

   
     // Read Variables
    prop.readValues( file );

   
    // Close File
    file.close();
  
    return prop.isEverythingRead();
    
}