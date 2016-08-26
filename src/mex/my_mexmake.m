% mex CXXFLAGS="\$CXXFLAGS -msse4.1 -falign-loops=16 -Wall -fopenmp" ...
%     DYLDFLAGS="\$DYLDFLAGS -fopenmp" mex_do_force_quad.cpp  
% mex CXXFLAGS="\$CXXFLAGS -msse4.1 -falign-loops=16 -Wall -fopenmp" -lgomp mex_do_force_quad.cpp 

mex CXXFLAGS="\$CXXFLAGS -O3  -msse4.1 -funroll-all-loops -falign-loops=16 -Wall -fopenmp" -lgomp mex/mex_saraspecquad.cpp 
