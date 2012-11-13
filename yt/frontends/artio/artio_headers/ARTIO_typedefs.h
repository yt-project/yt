#include <string>
#include <vector>
#include <cstdio>
#include <iostream>
#include <fstream>

#include "ARTIO_amr_data.hh"
#include "ARTIO_hydro_data.hh"

//... define the ARTIO base cell type to be of cell_locally_essential type
//... - this type allows moving between refinement levels
namespace ARTIO {
  namespace AMR{
    typedef ARTIO::AMR::cell_locally_essential<unsigned, float> ARTIO_cell;
    typedef ARTIO::AMR::level<ARTIO_cell> ARTIO_level;
    typedef ARTIO::AMR::tree< ARTIO_cell, ARTIO::AMR::level< ARTIO_cell > > ARTIO_tree;
  }
  namespace HYDRO {
    typedef ARTIO::HYDRO::data< ARTIO::AMR::ARTIO_tree > ARTIO_hydro_data;
  }
}
