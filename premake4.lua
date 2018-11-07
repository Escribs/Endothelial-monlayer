-- A solution contains projects, and defines the available configurations
solution "main"
   configurations { "Debug", "Release", "Cluster" }
 
   -- A project defines one build target
   project "Tester"
	kind "ConsoleApp"
	language "C++"
	files { "**.hpp", "**.cpp" }
	links {"gomp"} 
	buildoptions {"-std=c++11",
		    "-fopenmp"
			 
 				}
 
      configuration "Debug"
         defines { "DEBUG" }
         flags { "Symbols" }
	 buildoptions { "-isystem /usr/include/eigen3" } 
 
      configuration "Release"
         defines { "NDEBUG" }
         flags { "Optimize" }
	 buildoptions { "-isystem /usr/include/eigen3" }  
   
      configuration "Cluster"
         defines { "NDEBUG" }
         flags { "Optimize" } 
	 buildoptions { "-isystem ~/.local/lib/Eigen" } 
