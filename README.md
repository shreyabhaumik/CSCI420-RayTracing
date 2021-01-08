# **USC CSCI 420: Computer Graphics**  
#### \- taught by Dr. Jernej Barbič  

## **Programming Assignment 3: Ray Tracing**  

    Operating System: macOS 10.15
    Source Code Editor: Sublime Text, Version 3.2.2, Build 3211
    Programming Language: C++
    API: OpenGL (Core Profile)

### **ASSIGNMENT DETAILS:**
This assignment is to build a Ray Tracer which can handle opaque surfaces with lighting and shadows.  

- Look into the file ```Assignment-3-details.pdf```  
                OR
- Go to [![http://barbic.usc.edu/cs420-s20/assign3/index.html](http://barbic.usc.edu/cs420-s20/assign3/index.html)](http://barbic.usc.edu/cs420-s20/assign3/index.html)

### **HOW TO EXECUTE THE CODE (For macOS):**
1. Go to folder ```hw3-starterCode```
2. Compile using the command: ```make```.  
Note: To delete all the object files and executables:```make clean```
3. To execute type:  
    ```
    ./hw3 <input scenefile> [y(Y)/n(N) - to enable/disable SSAA Anti-aliasing] [y(Y)/n(N) - to enable/disable Soft Shadows] [numberOfReflectionsInInteger] [output jpegname]
    ```  
    Samples:  
    ```
    ./hw3 spheres.scene y y 1 outputspheres.jpg
	./hw3 table.scene y n 0 outputtable.jpg
	./hw3 test2.scene y
	./hw3 siggraph.scene
    ```  
  
If you want to only use antialiasing, and don't want to produce a jpg, you can leave out specifying the rest - it automatically takes all the rest as _no_.
After the ray tracing is complete, just press ESC to quit.  

### **MANDATORY FEATURES:**
1. Ray tracing triangles
2. Ray tracing spheres
3. Triangle Phong Shading
4. Sphere Phong Shading
5. Shadows rays
6. Still images

### **ADDITIONAL FEATURES(EXTRA CREDIT CONSIDERATIONS):**

1. **Recursive reflection**  
    Can produce as many reflection rays recursively. But, effect shows less for the higher recursions. (just pass the [numberOfReflectionsInInteger] as 3 or 4 or 5 or 0(for no reflection - only primary ray))

2. **Good antialiasing**  
    Used Supersampling Anti-aliasing method. Divided the pixel into 4 equal parts and sent rays through the center of each of those parts. Then averaged out the color contributed by the 4 rays by dividing the sum of the colors by 4.

3. **Soft shadows**  
    Created multiple point lights in a sphere centered on per point-light(provided with the starter code). Gave each subpoint-light their color value by dividing the color value of each point-lights by number of subpoint-lights created out of it. Then computed color of pixel produced by effect from all those sublights.


### **EXECUTION DURATION(FOR MANDATORY FEATURES):**
The core part for the respective scenes take the following durations:  
- test1.scene - 11 secs  
- test2.scene - 11 secs  
- spheres.scene - 11 secs  
- table.scene - 11 secs  
- siggraph.scene - 41 secs  

### **FILES AND FOLDERS LOCATIONS:**
- The ```.scene``` files are in the folder "hw3-starterCode".

### **EXHIBIT:**
![Reflection](Still-Images/000.jpg)  
REFLECTION(1 iteration) - test2.scene - 11 secs 	[```./hw3 test2.scene n n 1 000.jpg``` or ```./hw3 test2.scene n n 1```]  

![Antialiasing](Still-Images/001.jpg)  
ANTIALIASING - test2.scene - 11 secs 	[```./hw3 test2.scene y n 0 001.jpg``` or ```./hw3 test2.scene y```]  

![Soft Shadows](Still-Images/002.jpg)  
SOFTSHADOWS - table.scene - 42 secs 	[```./hw3 table.scene n y 0 002.jpg``` or ```./hw3 table.scene n y```]  

![1Reflection&Antialiasing](Still-Images/003.jpg)  
REFLECTION(1 iteration) + ANTIALIASING - spheres.scene - 11 secs 	[```./hw3 spheres.scene y n 1 003.jpg``` or ```./hw3 spheres.scene y n 1```]  

![3Reflections&Antialiasing](Still-Images/004.jpg)  
REFLECTION(3 iterations) + ANTIALIASING - spheres.scene - 11 secs 	[```./hw3 spheres.scene y n 3 004.jpg``` or ```./hw3 spheres.scene y n 3```]  

![1Reflection&SoftShadows](Still-Images/005.jpg)  
REFLECTION(1 iteration) + SOFTSHADOWS - spheres.scene - 11 secs 	[```./hw3 spheres.scene n y 1 005.jpg``` or ```./hw3 spheres.scene n y 1```]  

![Antialiasing&SoftShadows](Still-Images/006.jpg)  
ANTIALIASING + SOFTSHADOWS - table.scene - 2 mins 44 secs 	[```./hw3 table.scene y y 0 006.jpg``` or ```./hw3 table.scene y y```]  

![3Reflections&Antialiasing&SoftShadows-Spheres](Still-Images/007.jpg)  
REFLECTION(3 iterations) + ANTIALIASING + SOFTSHADOWS - spheres.scene - 14 secs 	[```./hw3 spheres.scene y y 3 007.jpg``` or ```./hw3 spheres.scene y y 3```]  

![3Reflections&Antialiasing&SoftShadows-Table](Still-Images/008.jpg)  
REFLECTION(3 iterations) + ANTIALIASING + SOFTSHADOWS - table.scene - 8 mins 57 secs 	[```./hw3 table.scene y y 3 008.jpg``` or ```./hw3 table.scene y y 3```]  

![3Reflections&Antialiasing&SoftShadows-Siggraph](Still-Images/009.jpg)  
REFLECTION(3 iterations) + ANTIALIASING + SOFTSHADOWS - siggraph.scene - 53 mins 9 secs 	[```./hw3 siggraph.scene y y 3 009.jpg``` or ```./hw3 siggraph.scene y y 3```]  