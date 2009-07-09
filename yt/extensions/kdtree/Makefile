# if you want to build this statically, you need to include fKD.f90 to the
# compile line, so pick the first of the two below. Otherwise, eliminate it, 
# like the second, for a shared object.


fKD: fKD.f90 fKD.v fKD_source.f90
#	Forthon --compile_first fKD_source --no2underscores --with-numpy -g fKD fKD.f90 fKD_source.f90
	Forthon -F gfortran --compile_first fKD_source --no2underscores --with-numpy --fopt "-O3" fKD fKD_source.f90

clean:
	rm -rf build fKDpy.a fKDpy.so
