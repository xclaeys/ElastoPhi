# ElastoPhi
Files of the CEMRACS project ElastoPhi


To compile : 
- create a folder "bin" or "build"
- In this folder do : cmake ..
- Do : make
- If you want to run tests, do : make test
- if you want to generate documentation, do : make doc
- the executable will be bin/src/main or build/src/main


To add a test :
- add the file in test/
- modify CMakeLists.txt and test/CMakeLists.txt 



To write latex in documentation :
- in-text forumulas : \f$ ... \f$
- unumbered centered displayed formulas : \f[ ... \f]
- special env : \f{env}{ ... \f}
 
To read the documentation :
- HTML : open index.html
- latex, do : make in bin/doc/latex and open bin/doc/latex/refman.pdf


