from grid_data_format import AthenaDistributedConverter
import sys
# Assumes that last input is the basename for the athena dataset.
# i.e. kh_3d_mhd_hlld_128_beta5000_sub_tanhd.0030
basename = sys.argv[-1]
converter = AthenaDistributedConverter(basename)
converter.convert()

# If you have units information, set up a conversion dictionary for
# each field.  Each key is the name of the field that Athena uses.
# Each value is what you have to multiply the raw output from Athena
# by to get cgs units.

# code_to_cgs = {'density':1.0e3,
# 	       'total_energy':1.0e-3,
# 	       'velocity_x':1.2345,
# 	       'velocity_y':1.2345,
# 	       'velocity_z':1.2345}

# converter = AthenaDistributedConverter(basename, field_conversions=code_to_cgs)
# converter.convert()

