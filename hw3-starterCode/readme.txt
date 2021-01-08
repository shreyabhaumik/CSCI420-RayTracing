Assignment #3: Ray tracing

FULL NAME: Shreya Bhaumik

Note:
used the starter code version at: http://barbic.usc.edu/cs420-s20/assign3/hw3-starterCode-vs2017.zip
got the error: error: use of undeclared identifier 'nullptr'
to resolve this, had to include -std=c++11 in the list of CXXFLAGS in the Makefile.

How to run the code:
After compiling using 'make',
execute: ./hw3 <input scenefile> [y(Y)/n(N) - to enable/disable SSAA Anti-aliasing] [y(Y)/n(N) - to enable/disable Soft Shadows] [numberOfReflectionsInInteger] [output jpegname]
example: ./hw3 spheres.scene y y 1 outputspheres.jpg
		 ./hw3 table.scene y n 0 outputtable.jpg
		 ./hw3 test2.scene y
		 ./hw3 siggraph.scene
If you want to only use antialiasing, and don't want to produce a jpg, you can leave out specifying the rest - it automatically takes all the rest as no.
After the ray tracing is complete, just press ESC to quit.

MANDATORY FEATURES
------------------

<Under "Status" please indicate whether it has been implemented and is
functioning correctly.  If not, please explain the current status.>

Feature:                                 Status: finish? (yes/no)
-------------------------------------    -------------------------
1) Ray tracing triangles                  			Yes

2) Ray tracing sphere                     			Yes

3) Triangle Phong Shading                 			Yes

4) Sphere Phong Shading                   			Yes

5) Shadows rays                           			Yes

6) Still images                           			Yes
   
7) Extra Credit (up to 20 points)
   1. Recursive reflection (10 points):
      Can produce as many reflection rays recursively. But, effect shows less for the higher recursions. (just pass the [numberOfReflectionsInInteger] as 3 or 4 or 5 or 0(for no reflection - only primary ray))
   2. Good antialiasing (10 points):
   	  Used Supersampling Anti-aliasing method. Divided the pixel into 4 equal parts and sent rays through the center of each of those parts. Then averaged out the color contributed by the 4 rays by dividing the sum of the colors by 4.
   3. Soft shadows (10 points):
   	  Created multiple point lights in a sphere centered on per point-light(provided with the starter code). Gave each subpoint-light their color value by dividing the color value of each point-lights by number of subpoint-lights created out of it. Then computed color of pixel produced by effect from all those sublights.


Still images:
000.jpg - REFLECTION(1 iteration) - test2.scene - 11 secs 	['./hw3 test2.scene n n 1 000.jpg' or './hw3 test2.scene n n 1']
001.jpg - ANTIALIASING - test2.scene - 11 secs 	['./hw3 test2.scene y n 0 001.jpg' or './hw3 test2.scene y']
002.jpg - SOFTSHADOWS - table.scene - 42 secs 	['./hw3 table.scene n y 0 002.jpg' or './hw3 table.scene n y']
003.jpg - REFLECTION(1 iteration) + ANTIALIASING - spheres.scene - 11 secs 	['./hw3 spheres.scene y n 1 003.jpg' or './hw3 spheres.scene y n 1']
004.jpg - REFLECTION(3 iteration) + ANTIALIASING - spheres.scene - 11 secs 	['./hw3 spheres.scene y n 3 004.jpg' or './hw3 spheres.scene y n 3']
005.jpg - REFLECTION(1 iteration) + SOFTSHADOWS - spheres.scene - 11 secs 	['./hw3 spheres.scene n y 1 005.jpg' or './hw3 spheres.scene n y 1']
006.jpg - ANTIALIASING + SOFTSHADOWS - table.scene - 2 mins 44 secs 	['./hw3 table.scene y y 0 006.jpg' or './hw3 table.scene y y']
007.jpg - REFLECTION(3 iteration) + ANTIALIASING + SOFTSHADOWS - spheres.scene - 14 secs 	['./hw3 spheres.scene y y 3 007.jpg' or './hw3 spheres.scene y y 3']
008.jpg - REFLECTION(3 iteration) + ANTIALIASING + SOFTSHADOWS - table.scene - 8 mins 57 secs 	['./hw3 table.scene y y 3 008.jpg' or './hw3 table.scene y y 3']
009.jpg - REFLECTION(3 iteration) + ANTIALIASING + SOFTSHADOWS - siggraph.scene - 53 mins 9 secs 	['./hw3 siggraph.scene y y 3 009.jpg' or './hw3 siggraph.scene y y 3']
Note: All the soft shadows are produced using 20 subpoint-lights per point-light.

Execution duration(for core part):
The core part for the respective scenes take the following durations:
test1.scene - 11 secs
test2.scene - 11 secs
spheres.scene - 11 secs
table.scene - 11 secs
siggraph.scene - 41 secs

Code edited/created by me are in two files:
hw3.cpp
myStructures.h

Thank you!