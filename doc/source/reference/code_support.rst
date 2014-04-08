
.. _code-support:

Code Support
============

Levels of Support for Various Codes
-----------------------------------

yt provides frontends to support several different simulation code formats 
as inputs.  Below is a list showing what level of support is provided for
each code.

|

+----------------------+------+-------+-------+------+---------+--------+-------+--------+
| Capability           | Enzo | Orion | FLASH | Nyx  | Piernik | Athena | Pluto | Chombo |
+======================+======+=======+=======+======+=========+========+=======+========+
| Fluid Quantities     |   Y  |   Y   |   Y   |  Y   |    Y    |   Y    |   Y   |    Y   |
+----------------------+------+-------+-------+------+---------+--------+-------+--------+
| Particles            |   Y  |   Y   |   Y   |  Y   |   N/A   |   N    |   N   |    N   |
+----------------------+------+-------+-------+------+---------+--------+-------+--------+
| Parameters           |   Y  |   Y   |   Y   |  Y   |    Y    |   Y    |   Y   |    Y   |
+----------------------+------+-------+-------+------+---------+--------+-------+--------+
| Units                |   Y  |   Y   |   Y   |  Y   |    Y    |   Y    |   Y   |    Y   |
+----------------------+------+-------+-------+------+---------+--------+-------+--------+
| Read on Demand       |   Y  |   Y   |   Y   |  Y   |    Y    |   Y    |   Y   |    Y   |
+----------------------+------+-------+-------+------+---------+--------+-------+--------+
| Load Raw Data        |   Y  |   Y   |   Y   |  Y   |    Y    |   Y    |   Y   |    Y   |
+----------------------+------+-------+-------+------+---------+--------+-------+--------+
| Part of test suite   |   Y  |   Y   |   Y   |  Y   |    N    |   N    |   N   |    Y   |
+----------------------+------+-------+-------+------+---------+--------+-------+--------+
| Level of Support     | Full | Full  | Full  | Full |  Full   |  Full  | Part  |  Part  |
+----------------------+------+-------+-------+------+---------+--------+-------+--------+

|

If you have a dataset from a code not yet supported, you can either 
input your data following :ref:`loading-numpy-array`, or help us by 
:ref:`creating_frontend` for this new format.

.. note::
   
   Support for additional codes, including particle codes (Gadget, Gasoline),
   octree AMR codes (ART, RAMSES), and additional patch AMR codes (Castro,
   Maestro, and other Boxlib codes) will be available in yt 3.0.  Currently yt
   3.0 is under active development, please stop by the mailing list if you want
   to get involved.

Future Codes to Support
-----------------------

A major overhaul of the code was required in order to cleanly support 
additional codes.  Development in the yt 3.x branch has begun and provides 
support for codes like: 
`RAMSES <http://irfu.cea.fr/Phocea/Vie_des_labos/Ast/ast_sstechnique.php?id_ast=904>`_, 
`ART (NMSU) <http://adsabs.harvard.edu/abs/1997ApJS..111...73K>`_, and 
`Gadget <http://www.mpa-garching.mpg.de/gadget/>`_.  Please switch to that 
version of yt for the most up-to-date support for those codes.

Additionally, in yt 3.0 the Boxlib formats have been unified and streamlined.
