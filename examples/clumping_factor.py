# Contributed by: Britton Smith (brittons@origins.colorado.edu)
#
# Calculates clumping factors over the entire simulation box.
# Arguments:
#     dataset_object: EnzoStaticOutput object.
#     fields: array containing the fields for the clumping factor.
#     exponent: array same size as fields of exponents for each of the fields.
#     weight: weighting field.
#
# Example: The standard density clumping factor (<rho^2>/<rho>^2) with 
#          mass-weighted averages.
#          clumping_factor(my_object,["Density"],[2],[CellMassMsun])
#          This is the same as:
#          clumping_factor(my_object,["Density","Density"],[1,1],[CellMassMsun])

import yt.lagos as lagos
from yt.funcs import get_pbar
import numpy as na

def clumping_factor(dataset_object,fields,exponent,weight):
    "Calculates clumping factor of array of fields with a weight."

    string_top = "<"
    string_bottom = "("
    for q in range(len(fields)):
        string_top = "%s %s^%s" % (string_top,fields[q],exponent[q])
        string_bottom = "%s <%s>^%s" % (string_bottom,fields[q],exponent[q])
    string_top = "%s >" % string_top
    string_bottom = "%s )" % string_bottom
    print "Clumping factor is %s / %s, weighted by %s." % (string_top,string_bottom,weight)

    product = 0.0
    weight_value = 0.0
    individual = [0.0 for q in fields]
    individual = na.array(individual)
    exponent = na.array(exponent)

    pb = get_pbar("Calculating clumping factor ", len(dataset_object.h.grids)) # Progress bar
    for grid in dataset_object.h.grids:
        pb.update(grid.id)
        this_product = grid[weight] * grid.child_mask
        for q in range(len(fields)):
            individual[q] += (grid[fields[q]] * grid[weight] * grid.child_mask).sum()
            this_product *= (grid[fields[q]]**exponent[q])
        product += this_product.sum()
        weight_value += (grid[weight] * grid.child_mask).sum()
        grid.clear_data()
    pb.finish()

    top = product * (weight_value**(exponent.sum()-1))
    bottom = (individual**exponent).prod()

    clumping_factor = top / bottom

    print "Clumping factor is %e." % clumping_factor
    return clumping_factor
