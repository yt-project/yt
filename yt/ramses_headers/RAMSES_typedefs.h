#include <string>
#include <vector>
#include <cstdio>
#include <iostream>
#include <fstream>

#include "RAMSES_amr_data.hh"
#include "RAMSES_hydro_data.hh"

//... define the RAMSES base cell type to be of cell_locally_essential type
//... - this type allows moving between refinement levels
namespace RAMSES {
  namespace AMR{
    typedef RAMSES::AMR::cell_locally_essential<unsigned, float> RAMSES_cell;
    typedef RAMSES::AMR::level<RAMSES_cell> RAMSES_level;
    typedef RAMSES::AMR::tree< RAMSES_cell, RAMSES::AMR::level< RAMSES_cell > > RAMSES_tree;
  }
  namespace HYDRO {
    typedef RAMSES::HYDRO::data< RAMSES::AMR::RAMSES_tree > RAMSES_hydro_data;
  }
}
