# ElastoPhi
Files of the CEMRACS project ElastoPhi

## Compilation

- create a folder "bin" or "build"
- In this folder do : cmake ..
- Do : make   (add -jN to use N proc)
- If you want to run tests, do : make test
- if you want to generate documentation, do : make doc
- the executable will be bin/src/main or build/src/main

Note : -DCMAKE_BUILD_TYPE=debug or release

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
To generate the plots of the multi-compress results, in postprocessing directory do:
python graphes_output_compression.py ../output/output_compression_matrice....txt 

## Testing

### To add a test
- add the file in test/
- modify CMakeLists.txt and test/CMakeLists.txt


## Documentation

### To write documentation :


###To read the documentation :
- HTML : open index.html
- latex, do : make in bin/doc/latex and open bin/doc/latex/refman.pdf

###To write latex in documentation :
- in-text forumulas : \f$ ... \f$
- unumbered centered displayed formulas : \f[ ... \f]
- special env : \f{env}{ ... \f}


