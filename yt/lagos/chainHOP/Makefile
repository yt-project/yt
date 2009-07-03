# if you want to build this statically, you need to include fKD.f90 to the
# compile line, so pick the first of the two below. Otherwise, eliminate it, 
# like the second, for a shared object.


fchainHOP: fchainHOP.f90 fchainHOP.v fKD_source.f90
	Forthon --compile_first fKD_source --no2underscores --with-numpy -g fchainHOP fKD_source.f90

clean:
	rm -rf build fchainHOPpy.a fchainHOPpy.so
