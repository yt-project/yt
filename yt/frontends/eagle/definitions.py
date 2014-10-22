"""
EAGLE definitions




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

eaglenetwork_ions = \
    ('electron', 'H1', 'H2', 'H_m', 'He1', 'He2','He3', 'C1',\
     'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C_m', 'N1', 'N2', \
     'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'O1', 'O2', 'O3', \
     'O4', 'O5', 'O6', 'O7', 'O8', 'O9', 'O_m', 'Ne1', 'Ne2',\
     'Ne3', 'Ne4', 'Ne5', 'Ne6', 'Ne7', 'Ne8', 'Ne9', 'Ne10',\
     'Ne11', 'Mg1', 'Mg2', 'Mg3', 'Mg4', 'Mg5', 'Mg6', 'Mg7',\
     'Mg8', 'Mg9', 'Mg10', 'Mg11', 'Mg12', 'Mg13', 'Si1', 'Si2',\
     'Si3', 'Si4', 'Si5', 'Si6', 'Si7', 'Si8', 'Si9', 'Si10',\
     'Si11', 'Si12', 'Si13', 'Si14', 'Si15', 'Si16', 'Si17',\
     'Ca1', 'Ca2', 'Ca3', 'Ca4', 'Ca5', 'Ca6', 'Ca7', 'Ca8',\
     'Ca9', 'Ca10', 'Ca11', 'Ca12', 'Ca13', 'Ca14', 'Ca15',\
     'Ca16', 'Ca17', 'Ca18', 'Ca19', 'Ca20', 'Ca21', 'Fe1',\
     'Fe2', 'Fe3', 'Fe4', 'Fe5', 'Fe6', 'Fe7', 'Fe8', 'Fe9',\
     'Fe10', 'Fe11', 'Fe12', 'Fe13', 'Fe14', 'Fe15', 'Fe16',\
     'Fe17', 'Fe18', 'Fe19', 'Fe20', 'Fe21', 'Fe22', 'Fe23',\
     'Fe24', 'Fe25', 'Fe25', 'Fe27',)

eaglenetwork_ion_lookup = {ion:index for index, ion in enumerate(eaglenetwork_ions)}
