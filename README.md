# ElastoPhi
Files of the CEMRACS project ElastoPhi

## Input files
Some executables require an input file with:  
Eta : ...  
Epsilon : ...  
Data_path : ...  
Output_path : ...  
Mesh_name : ...  
Matrix_name : ...  
Some executables don't necessarily require all these parameters.

Example of input file:  
Eta : 1  
Epsilon : 0.9  
Data_path : data  
Output_path : output  
Mesh_name : maillage1994Fracs.txt  
Matrix_name : matrice1994Fracs.txt  

## Postprocessing
To generate plots, python scripts usually just need the data file in input

## Examples of results

![Cluster tree for fractures](maillage450Fracs.gif?raw=true)  ![Cluster tree for faults](maillage5364FracsTriangles.gif?raw=true)

![Compression](compression.gif?raw=true)
