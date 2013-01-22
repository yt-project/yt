# if you want to build this statically, you need to include fKD.f90 to the
# compile line, so pick the first of the two below. Otherwise, eliminate it, 
# like the second, for a shared object.

ifdef FORTHON_EXE
FORTHON=$(FORTHON_EXE)
else
FORTHON=Forthon
endif

fKD: fKD.f90 fKD.v fKD_source.f90
#	Forthon --compile_first fKD_source --no2underscores -g fKD fKD.f90 fKD_source.f90
	@echo "Using $(FORTHON) ($(FORTHON_EXE))"
	$(FORTHON) -F gfortran --compile_first fKD_source --no2underscores --fopt "-O3" fKD fKD_source.f90
	mv build/lib*/fKDpy.so .

clean:
	rm -rf build fKDpy.a fKDpy.so
