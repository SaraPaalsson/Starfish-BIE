include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)


## MEX functions
matlab_add_mex(
	NAME mex_saraspecquad
	SRC mex_saraspecquad.cpp
)
