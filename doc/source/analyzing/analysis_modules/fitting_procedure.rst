.. _fitting_procedure:

Procedure for Generating Fits
=============================
.. sectionauthor:: Hilary Egan <hilary.egan@colorado.edu>

To generate a fit for a spectrum :py:func:`generate_total_fit()` is called.
This function controls the identification of line complexes, the fit
of a series of absorption lines for each appropriate species, checks of
those fits, and returns the results of the fits.


Finding Line Complexes
----------------------
Line complexes are found using the :py:func:`find_complexes` function. The
process by which line complexes are found involves walking through
the array of flux in order from minimum to maximum wavelength, and finding
series of spatially contiguous cells whose flux is less than some limit.
These regions are then checked in terms of an additional flux limit and size.
The bounds of all the passing regions are then listed and returned. Those
bounds that cover an exceptionally large region of wavelength space will be
broken up if a suitable cut point is found. This method is only appropriate
for noiseless spectra.

The optional parameter **complexLim** (default = 0.999), controls the limit
that triggers the identification of a spatially contiguous region of flux
that could be a line complex. This number should be very close to 1 but not
exactly equal. It should also be at least an order of magnitude closer to 1
than the later discussed **fitLim** parameter, because a line complex where
the flux of the trough is very close to the flux of the edge can be incredibly
unstable when optimizing.

The **fitLim** parameter controls what is the maximum flux that the trough
of the region can have and still be considered a line complex. This 
effectively controls the sensitivity to very low column absorbers. Default
value is **fitLim** = 0.99. If a region is identified where the flux of the trough
is greater than this value, the region is simply ignored.

The **minLength** parameter controls the minimum number of array elements 
that an identified region must have. This value must be greater than or
equal to 3 as there are a minimum of 3 free parameters that must be fit.
Default is **minLength** = 3.

The **maxLength** parameter controls the maximum number of array elements
that an identified region can have before it is split into separate regions.
Default is **maxLength** = 1000. This should be adjusted based on the 
resolution of the spectrum to remain appropriate. The value correspond
to a wavelength of roughly 50 angstroms. 

The **splitLim** parameter controls how exceptionally large regions are split.
When such a region is identified by having more array elements than
**maxLength**, the point of maximum flux (or minimum absorption) in the 
middle two quartiles is identified. If that point has a flux greater than
or equal to **splitLim**, then two separate complexes are created: one from
the lower wavelength edge to the minimum absorption point and the other from
the minimum absorption point to the higher wavelength edge. The default
value is **splitLim** =.99, but it should not drastically affect results, so
long as the value is reasonably close to 1.


Fitting a Line Complex
----------------------

After a complex is identified, it is fitted by iteratively adding and 
optimizing a set of Voigt Profiles for a particular species until the
region is considered successfully fit. The optimizing is accomplished
using scipy's least squares optimizer. This requires an initial estimate
of the parameters to be fit (column density, b-value, redshift) for each
line.

Each time a line is added, the guess of the parameters is based on
the difference between the line complex and the fit so far. For the first line
this just means the initial guess is based solely on the flux of the line
complex. The column density is given by the initial column density given
in the species parameters dictionary. If the line is saturated (some portion
of the flux with a value less than .1) than the larger initial column density
guess is chosen. If the flux is relatively high (all values >.9) than the
smaller initial guess is given. These values are chosen to make optimization
faster and more stable by being closer to the actual value, but the final
results of fitting should not depend on them as they merely provide a
starting point. 

After the parameters for a line are optimized for the first time, the 
optimized parameters are then used for the initial guess on subsequent 
iterations with more lines. 

The complex is considered successfully fit when the sum of the squares of 
the difference between the flux generated from the fit and the desired flux
profile is less than **errBound**. **errBound** is related to the optional
parameter to :py:func:`generate_total_fit()`, **maxAvgError** by the number
of array elements in the region such that **errBound** = number of elements *
**maxAvgError**.

There are several other conditions under which the cycle of adding and 
optimizing lines will halt. If the error of the optimized fit from adding
a line is an order of magnitude worse than the error of the fit without
that line, then it is assumed that the fitting has become unstable and 
the latest line is removed. Lines are also prevented from being added if
the total number of lines is greater than the number of elements in the flux
array being fit divided by 3. This is because there must not be more free
parameters in a fit than the number of points to constrain them. 


Checking Fit Results
--------------------

After an acceptable fit for a region is determined, there are several steps
the algorithm must go through to validate the fits. 

First, the parameters must be in a reasonable range. This is a check to make 
sure that the optimization did not become unstable and generate a fit that
diverges wildly outside the region where the fit was performed. This way, even
if particular complex cannot be fit, the rest of the spectrum fitting still
behaves as expected. The range of acceptability for each parameter is given
in the species parameter dictionary. These are merely broad limits that will
prevent numerical instability rather than physical limits.

In cases where a single species generates multiple lines (as in the OVI 
doublet), the fits are then checked for higher wavelength lines. Originally
the fits are generated only considering the lowest wavelength fit to a region.
This is because we perform the fitting of complexes in order from the lowest
wavelength to the highest, so any contribution to a complex being fit must
come from the lower wavelength as the higher wavelength contributions would
already have been subtracted out after fitting the lower wavelength. 

Saturated Lyman Alpha Fitting Tools
-----------------------------------

In cases where a large or saturated line (there exists a point in the complex
where the flux is less than .1) fails to be fit properly at first pass, a
more robust set of fitting tools is used to try and remedy the situation.
The basic approach is to simply try a much wider range of initial parameter
guesses in order to find the true optimization minimum, rather than getting
stuck in a local minimum. A set of hard coded initial parameter guesses
for Lyman alpha lines is given by the function :py:func:`get_test_lines`. 
Also included in these parameter guesses is an an initial guess of a high
column cool line overlapping a lower column warm line, indictive of a 
broad Lyman alpha (BLA) absorber.
