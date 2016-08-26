% mex CXXFLAGS="\$CXXFLAGS -msse4.1 -falign-loops=16 -Wall -fopenmp" ...
%     DYLDFLAGS="\$DYLDFLAGS -fopenmp" mex_do_force_quad.cpp  
% mex CXXFLAGS="\$CXXFLAGS -msse4.1 -falign-loops=16 -Wall -fopenmp" -lgomp mex_do_force_quad.cpp 

mex CXXFLAGS="\$CXXFLAGS -O3  -msse4.1 -Wall " mex_saraspecquad.cpp 
