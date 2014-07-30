
.. _code-support:

Code Support
============

Levels of Support for Various Codes
-----------------------------------

yt provides frontends to support several different simulation code formats 
as inputs.  Below is a list showing what level of support is provided for
each code.

|

+-----------------------+------------+-----------+------------+-------+---------+----------+------------+----------+ 
| Code ▼ \ Capability ► |   Fluid    | Particles | Parameters | Units | Read on | Load Raw |  Part of   | Level of | 
|                       | Quantities |           |            |       | Demand  |   Data   | test suite | Support  |
+-----------------------+------------+-----------+------------+-------+---------+----------+------------+----------+ 
| Athena                |     Y      |     N     |      Y     |   Y   |    Y    |    Y     |     N      |   Full   |
+-----------------------+------------+-----------+------------+-------+---------+----------+------------+----------+ 
| Castro                |     Y      |     Y     |   Partial  |   Y   |    Y    |    Y     |     N      |   Full   |
+-----------------------+------------+-----------+------------+-------+---------+----------+------------+----------+ 
| Chombo                |     Y      |     N     |      Y     |   Y   |    Y    |    Y     |     Y      | Partial  |
+-----------------------+------------+-----------+------------+-------+---------+----------+------------+----------+ 
| Enzo                  |     Y      |     Y     |      Y     |   Y   |    Y    |    Y     |     Y      |   Full   |
+-----------------------+------------+-----------+------------+-------+---------+----------+------------+----------+ 
| FLASH                 |     Y      |     Y     |      Y     |   Y   |    Y    |    Y     |     Y      |   Full   |
+-----------------------+------------+-----------+------------+-------+---------+----------+------------+----------+ 
| Maestro               |   Y [#f1_] |     N     |      Y     |   Y   |    Y    |    Y     |     N      | Partial  |
+-----------------------+------------+-----------+------------+-------+---------+----------+------------+----------+ 
| Nyx                   |     Y      |     Y     |      Y     |   Y   |    Y    |    Y     |     Y      |   Full   |
+-----------------------+------------+-----------+------------+-------+---------+----------+------------+----------+ 
| Orion                 |     Y      |     Y     |      Y     |   Y   |    Y    |    Y     |     Y      |   Full   |
+-----------------------+------------+-----------+------------+-------+---------+----------+------------+----------+ 
| Piernik               |     Y      |    N/A    |      Y     |   Y   |    Y    |    Y     |     N      |   Full   |
+-----------------------+------------+-----------+------------+-------+---------+----------+------------+----------+ 
| Pluto                 |     Y      |     N     |      Y     |   Y   |    Y    |    Y     |     N      | Partial  |
+-----------------------+------------+-----------+------------+-------+---------+----------+------------+----------+ 

.. [#f1] one-dimensional base-state not read in currently.

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
