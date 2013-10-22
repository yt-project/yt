from scipy import optimize
import numpy as na
import h5py
from yt.analysis_modules.absorption_spectrum.absorption_line \
        import voigt


def generate_total_fit(x, fluxData, orderFits, speciesDicts, 
        minError=1E-5, complexLim=.999,
        fitLim=.99, minLength=3, 
        maxLength=1000, splitLim=.99,
        output_file=None):

    """
    This function is designed to fit an absorption spectrum by breaking 
    the spectrum up into absorption complexes, and iteratively adding
    and optimizing voigt profiles to each complex.

    Parameters
    ----------
    x : (N) ndarray
        1d array of wavelengths
    fluxData : (N) ndarray
        array of flux corresponding to the wavelengths given
        in x. (needs to be the same size as x)
    orderFits : list
        list of the names of the species in the order that they 
        should be fit. Names should correspond to the names of the species
        given in speciesDicts. (ex: ['lya','OVI'])
    speciesDicts : dictionary
        Dictionary of dictionaries (I'm addicted to dictionaries, I
        confess). Top level keys should be the names of all the species given
        in orderFits. The entries should be dictionaries containing all 
        relevant parameters needed to create an absorption line of a given 
        species (f,Gamma,lambda0) as well as max and min values for parameters
        to be fit
    complexLim : float, optional
        Maximum flux to start the edge of an absorption complex. Different 
        from fitLim because it decides extent of a complex rather than 
        whether or not a complex is accepted. 
    fitLim : float,optional
        Maximum flux where the level of absorption will trigger 
        identification of the region as an absorption complex. Default = .98.
        (ex: for a minSize=.98, a region where all the flux is between 1.0 and
        .99 will not be separated out to be fit as an absorbing complex, but
        a region that contains a point where the flux is .97 will be fit
        as an absorbing complex.)
    minLength : int, optional
        number of cells required for a complex to be included. 
        default is 3 cells.
    maxLength : int, optional
        number of cells required for a complex to be split up. Default
        is 1000 cells.
    splitLim : float, optional
        if attempting to split a region for being larger than maxlength
        the point of the split must have a flux greater than splitLim 
        (ie: absorption greater than splitLim). Default= .99.
    output_file : string, optional
        location to save the results of the fit. 

    Returns
    -------
    allSpeciesLines : dictionary
        Dictionary of dictionaries representing the fit lines. 
        Top level keys are the species given in orderFits and the corresponding
        entries are dictionaries with the keys 'N','b','z', and 'group#'. 
        Each of these corresponds to a list of the parameters for every
        accepted fitted line. (ie: N[0],b[0],z[0] will create a line that
        fits some part of the absorption spectrum). 'group#' is a similar list
        but identifies which absorbing complex each line belongs to. Lines
        with the same group# were fit at the same time. group#'s do not
        correlate between species (ie: an lya line with group number 1 and
        an OVI line with group number 1 were not fit together and do
        not necessarily correspond to the same region)
    yFit : (N) ndarray
        array of flux corresponding to the combination of all fitted
        absorption profiles. Same size as x.
    """

    #Empty dictionary for fitted lines
    allSpeciesLines = {}

    #Wavelength of beginning of array, wavelength resolution
    x0,xRes=x[0],x[1]-x[0]

    #Empty fit without any lines
    yFit = na.ones(len(fluxData))

    #Force the first and last flux pixel to be 1 to prevent OOB
    fluxData[0]=1
    fluxData[-1]=1

    #Find all regions where lines/groups of lines are present
    cBounds = _find_complexes(x, fluxData, fitLim=fitLim,
            complexLim=complexLim, minLength=minLength,
            maxLength=maxLength, splitLim=splitLim)

    #Fit all species one at a time in given order from low to high wavelength
    for species in orderFits:
        speciesDict = speciesDicts[species]
        speciesLines = {'N':na.array([]),
                        'b':na.array([]),
                        'z':na.array([]),
                        'group#':na.array([])}

        #Set up wavelengths for species
        initWl = speciesDict['wavelength'][0]

        for b_i,b in enumerate(cBounds):
            xBounded=x[b[1]:b[2]]
            yDatBounded=fluxData[b[1]:b[2]]
            yFitBounded=yFit[b[1]:b[2]]

            #Find init redshift
            z=(xBounded[yDatBounded.argmin()]-initWl)/initWl

            #Check if any flux at partner sites
            if not _line_exists(speciesDict['wavelength'],
                    fluxData,z,x0,xRes,fitLim): 
                continue 

            #Fit Using complex tools
            newLinesP,flag=_complex_fit(xBounded,yDatBounded,yFitBounded,
                    z,fitLim,minError*(b[2]-b[1]),speciesDict)

            #Check existence of partner lines if applicable
            if len(speciesDict['wavelength']) != 1:
                newLinesP = _remove_unaccepted_partners(newLinesP, x, fluxData, 
                        b, minError*(b[2]-b[1]),
                        x0, xRes, speciesDict)

            #If flagged as a bad fit, species is lyman alpha,
            #   and it may be a saturated line, use special tools
            if flag and species=='lya' and min(yDatBounded)<.1:
                newLinesP=_large_flag_fit(xBounded,yDatBounded,
                        yFitBounded,z,speciesDict,
                        minSize,minError*(b[2]-b[1]))

            #Adjust total current fit
            yFit=yFit*_gen_flux_lines(x,newLinesP,speciesDict)

            #Add new group to all fitted lines
            if na.size(newLinesP)>0:
                speciesLines['N']=na.append(speciesLines['N'],newLinesP[:,0])
                speciesLines['b']=na.append(speciesLines['b'],newLinesP[:,1])
                speciesLines['z']=na.append(speciesLines['z'],newLinesP[:,2])
                groupNums = b_i*na.ones(na.size(newLinesP[:,0]))
                speciesLines['group#']=na.append(speciesLines['group#'],groupNums)

        allSpeciesLines[species]=speciesLines

    if output_file:
        _output_fit(allSpeciesLines, output_file)

    return (allSpeciesLines,yFit)

def _complex_fit(x, yDat, yFit, initz, minSize, errBound, speciesDict, 
        initP=None):
    """ Fit an absorption complex by iteratively adding and optimizing
    voigt profiles.
    
    A complex is defined as a region where some number of lines may be present,
    or a region of non zero of absorption. Lines are iteratively added
    and optimized until the difference between the flux generated using
    the optimized parameters has a least squares difference between the 
    desired flux profile less than the error bound.

    Parameters
    ----------
    x : (N) ndarray
        array of wavelength
    ydat : (N) ndarray
        array of desired flux profile to be fitted for the wavelength
        space given by x. Same size as x.
    yFit : (N) ndarray
        array of flux profile fitted for the wavelength
        space given by x already. Same size as x.
    initz : float
        redshift to try putting first line at 
        (maximum absorption for region)
    minsize : float
        minimum absorption allowed for a line to still count as a line
        given in normalized flux (ie: for minSize=.9, only lines with minimum
        flux less than .9 will be fitted)
    errbound : float
        maximum total error allowed for an acceptable fit
    speciesDict : dictionary
        dictionary containing all relevant parameters needed
        to create an absorption line of a given species (f,Gamma,lambda0)
        as well as max and min values for parameters to be fit
    initP : (,3,) ndarray
        initial guess to try for line parameters to fit the region. Used
        by large_flag_fit. Default = None, and initial guess generated
        automatically.

    Returns
    -------
    linesP : (3,) ndarray
        Array of best parameters if a good enough fit is found in 
        the form [[N1,b1,z1], [N2,b2,z2],...]
    flag : bool
        boolean value indicating the success of the fit (True if unsuccessful)
    """

    #Setup initial line guesses
    if initP==None: #Regular fit
        initP = [0,0,0] 
        if min(yDat)<.5: #Large lines get larger initial guess 
            initP[0] = 10**16
        elif min(yDat)>.9: #Small lines get smaller initial guess
            initP[0] = 10**12.5
        else:
            initP[0] = speciesDict['init_N']
        initP[1] = speciesDict['init_b']
        initP[2]=initz
        initP=na.array([initP])

    linesP = initP

    #For generating new z guesses
    wl0 = speciesDict['wavelength'][0]

    #Check if first line exists still
    if min(yDat-yFit+1)>minSize: 
        return [],False
    
    #Values to proceed through first run
    errSq,prevErrSq=1,1000

    while True:
        #Initial parameter guess from joining parameters from all lines
        #   in lines into a single array
        initP = linesP.flatten()

        #Optimize line
        fitP,success=optimize.leastsq(_voigt_error,initP,
                args=(x,yDat,yFit,speciesDict),
                epsfcn=1E-10,maxfev=1000)

        #Set results of optimization
        linesP = na.reshape(fitP,(-1,3))

        #Generate difference between current best fit and data
        yNewFit=_gen_flux_lines(x,linesP,speciesDict)
        dif = yFit*yNewFit-yDat

        #Sum to get idea of goodness of fit
        errSq=sum(dif**2)

        #If good enough, break
        if errSq < errBound: 
            break

        #If last fit was worse, reject the last line and revert to last fit
        if errSq > prevErrSq*10:
            #If its still pretty damn bad, cut losses and try flag fit tools
            if prevErrSq >1E2*errBound and speciesDict['name']=='HI lya':
                return [],True
            else:
                yNewFit=_gen_flux_lines(x,prevLinesP,speciesDict)
                break

        #If too many lines 
        if na.shape(linesP)[0]>8 or na.size(linesP)+3>=len(x):
            #If its fitable by flag tools and still bad, use flag tools
            if errSq >1E2*errBound and speciesDict['name']=='HI lya':
                return [],True
            else:
                break 

        #Store previous data in case reject next fit
        prevErrSq = errSq
        prevLinesP = linesP


        #Set up initial condition for new line
        newP = [0,0,0] 
        if min(dif)<.1:
            newP[0]=10**12
        elif min(dif)>.9:
            newP[0]=10**16
        else:
            newP[0]=10**14
        newP[1] = speciesDict['init_b']
        newP[2]=(x[dif.argmax()]-wl0)/wl0
        linesP=na.append(linesP,[newP],axis=0)


    #Check the parameters of all lines to see if they fall in an
    #   acceptable range, as given in dict ref
    remove=[]
    for i,p in enumerate(linesP):
        check=_check_params(na.array([p]),speciesDict)
        if check: 
            remove.append(i)
    linesP = na.delete(linesP,remove,axis=0)

    return linesP,False

def _large_flag_fit(x, yDat, yFit, initz, speciesDict, minSize, errBound):
    """
    Attempts to more robustly fit saturated lyman alpha regions that have
    not converged to satisfactory fits using the standard tools.

    Uses a preselected sample of a wide range of initial parameter guesses
    designed to fit saturated lines (see get_test_lines).

    Parameters
    ----------
    x : (N) ndarray
        array of wavelength
    ydat : (N) ndarray
        array of desired flux profile to be fitted for the wavelength
        space given by x. Same size as x.
    yFit : (N) ndarray
        array of flux profile fitted for the wavelength
        space given by x already. Same size as x.
    initz : float
        redshift to try putting first line at 
        (maximum absorption for region)
    speciesDict : dictionary
        dictionary containing all relevant parameters needed
        to create an absorption line of a given species (f,Gamma,lambda0)
        as well as max and min values for parameters to be fit
    minsize : float
        minimum absorption allowed for a line to still count as a line
        given in normalized flux (ie: for minSize=.9, only lines with minimum
        flux less than .9 will be fitted)
    errbound : float
        maximum total error allowed for an acceptable fit

    Returns
    -------
    bestP : (3,) ndarray
        array of best parameters if a good enough fit is found in 
        the form [[N1,b1,z1], [N2,b2,z2],...]  
    """

    #Set up some initial line guesses
    lineTests = _get_test_lines(initz)

    #Keep track of the lowest achieved error
    bestError = 1000 

    #Iterate through test line guesses
    for initLines in lineTests:
        if initLines[1,0]==0:
            initLines = na.delete(initLines,1,axis=0)

        #Do fitting with initLines as first guess
        linesP,flag=_complex_fit(x,yDat,yFit,initz,
                minSize,errBound,speciesDict,initP=initLines)

        #Find error of last fit
        yNewFit=_gen_flux_lines(x,linesP,speciesDict)
        dif = yFit*yNewFit-yDat
        errSq=sum(dif**2)

        #If error lower, keep track of the lines used to make that fit
        if errSq < bestError:
            bestError = errSq
            bestP = linesP

    if bestError>10*errBound*len(x): 
        return []
    else:
        return bestP

def _get_test_lines(initz):
    """
    Returns a 3d numpy array of lines to test as initial guesses for difficult
    to fit lyman alpha absorbers that are saturated. 
    
    The array is 3d because
    the first dimension gives separate initial guesses, the second dimension
    has multiple lines for the same guess (trying a broad line plus a 
    saturated line) and the 3d dimension contains the 3 fit parameters (N,b,z)

    Parameters
    ----------
    initz : float
        redshift to give all the test lines

    Returns
    -------
    testP : (,3,) ndarray
        numpy array of the form 
        [[[N1a,b1a,z1a], [N1b,b1b,z1b]], [[N2a,b2,z2a],...] ...]
    """

    #Set up a bunch of empty lines
    testP = na.zeros((10,2,3))

    testP[0,0,:]=[1E18,20,initz]
    testP[1,0,:]=[1E18,40,initz]
    testP[2,0,:]=[1E16,5, initz]
    testP[3,0,:]=[1E16,20,initz]
    testP[4,0,:]=[1E16,80,initz]

    testP[5,0,:]=[1E18,20,initz]
    testP[6,0,:]=[1E18,40,initz]
    testP[7,0,:]=[1E16,5, initz]
    testP[8,0,:]=[1E16,20,initz]
    testP[9,0,:]=[1E16,80,initz]

    testP[5,1,:]=[1E13,100,initz]
    testP[6,1,:]=[1E13,100,initz]
    testP[7,1,:]=[1E13,100,initz]
    testP[8,1,:]=[1E13,100,initz]
    testP[9,1,:]=[1E13,100,initz]

    return testP

def _get_bounds(z, b, wl, x0, xRes):
    """ 
    Gets the indices of range of wavelength that the wavelength wl is in 
    with the size of some initial wavelength range.

    Used for checking if species with multiple lines (as in the OVI doublet)
    fit all lines appropriately.

    Parameters
    ----------
    z : float
        redshift
    b : (3) ndarray/list
        initial bounds in form [i0,i1,i2] where i0 is the index of the 
        minimum flux for the complex, i1 is index of the lower wavelength 
        edge of the complex, and i2 is the index of the higher wavelength
        edge of the complex.
    wl : float
        unredshifted wavelength of the peak of the new region 
    x0 : float
        wavelength of the index 0
    xRes : float
        difference in wavelength for two consecutive indices
    
    Returns
    -------
    indices : (2) tuple
        Tuple (i1,i2) where i1 is the index of the lower wavelength bound of 
        the new region and i2 is the index of the higher wavelength bound of
        the new region
    """

    r=[-b[1]+100+b[0],b[2]+100-b[0]]
    redWl = (z+1)*wl
    iRedWl=int((redWl-x0)/xRes)
    indices = (iRedWl-r[0],iRedWl+r[1])

    return indices

def _remove_unaccepted_partners(linesP, x, y, b, errBound, 
        x0, xRes, speciesDict):
    """
    Given a set of parameters [N,b,z] that form multiple lines for a given
    species (as in the OVI doublet), remove any set of parameters where
    not all transition wavelengths have a line that matches the fit.

    (ex: if a fit is determined based on the first line of the OVI doublet,
    but the given parameters give a bad fit of the wavelength space of
    the second line then that set of parameters is removed from the array
    of line parameters.)

    Parameters
    ----------
    linesP : (3,) ndarray
        array giving sets of line parameters in 
        form [[N1, b1, z1], ...]
    x : (N) ndarray
        wavelength array [nm]
    y : (N) ndarray
        normalized flux array of original data
    b : (3) tuple/list/ndarray
        indices that give the bounds of the original region so that another 
        region of similar size can be used to determine the goodness
        of fit of the other wavelengths
    errBound : float
        size of the error that is appropriate for a given region, 
        adjusted to account for the size of the region.

    Returns
    -------
    linesP : (3,) ndarray
        array similar to linesP that only contains lines with
        appropriate fits of all transition wavelengths.
    """

    #List of lines to remove
    removeLines=[]

    #Iterate through all sets of line parameters
    for i,p in enumerate(linesP):

        #iterate over all transition wavelengths
        for wl in speciesDict['wavelength']:

            #Get the bounds of a similar sized region around the
            #   appropriate wavelength, and then get the appropriate
            #   region of wavelength and flux
            lb = _get_bounds(p[2],b,wl,x0,xRes)
            xb,yb=x[lb[0]:lb[1]],y[lb[0]:lb[1]]

            #Generate a fit and find the difference to data
            yFitb=_gen_flux_lines(xb,na.array([p]),speciesDict)
            dif =yb-yFitb

            #Only counts as an error if line is too big ---------------<
            dif = [k for k in dif if k>0]
            err = sum(dif)

            #If the fit is too bad then add the line to list of removed lines
            if err > errBound*1E2:
                removeLines.append(i)
                break

    #Remove all bad line fits
    linesP = na.delete(linesP,removeLines,axis=0)

    return linesP 



def _line_exists(wavelengths, y, z, x0, xRes,fluxMin):
    """For a group of lines finds if the there is some change in flux greater
    than some minimum at the same redshift with different initial wavelengths

    Parameters
    ----------
    wavelengths : (N) ndarray
        array of initial wavelengths to check
    y : (N) ndarray
        flux array to check
    x0 : float
        wavelength of the first value in y
    xRes : float
        difference in wavelength between consecutive cells in flux array
    fluxMin : float
        maximum flux to count as a line existing. 

    Returns
    -------

    flag : boolean 
        value indicating whether all lines exist. True if all lines exist
    """

    #Iterate through initial wavelengths
    for wl in wavelengths:
        #Redshifted wavelength
        redWl = (z+1)*wl

        #Index of the redshifted wavelength
        indexRedWl = (redWl-x0)/xRes

        #Check to see if even in flux range
        if indexRedWl > len(y):
            return False

        #Check if surpasses minimum absorption bound
        if y[int(indexRedWl)]>fluxMin:
            return False

    return True

def _find_complexes(x, yDat, complexLim=.999, fitLim=.99,
        minLength =3, maxLength=1000, splitLim=.99):
    """Breaks up the wavelength space into groups
    where there is some absorption. 

    Parameters
    ----------
    x : (N) ndarray
        array of wavelengths
    yDat : (N) ndarray
        array of flux corresponding to the wavelengths given
        in x. (needs to be the same size as x)
    complexLim : float, optional
        Maximum flux to start the edge of an absorption complex. Different 
        from fitLim because it decides extent of a complex rather than 
        whether or not a complex is accepted. 
    fitLim : float,optional
        Maximum flux where the level of absorption will trigger 
        identification of the region as an absorption complex. Default = .98.
        (ex: for a minSize=.98, a region where all the flux is between 1.0 and
        .99 will not be separated out to be fit as an absorbing complex, but
        a region that contains a point where the flux is .97 will be fit
        as an absorbing complex.)
    minLength : int, optional
        number of cells required for a complex to be included. 
        default is 3 cells.
    maxLength : int, optional
        number of cells required for a complex to be split up. Default
        is 1000 cells.
    splitLim : float, optional
        if attempting to split a region for being larger than maxlength
        the point of the split must have a flux greater than splitLim 
        (ie: absorption greater than splitLim). Default= .99.

    Returns
    -------
    cBounds : (3,) 
        list of bounds in the form [[i0,i1,i2],...] where i0 is the 
        index of the maximum flux for a complex, i1 is the index of the
        beginning of the complex, and i2 is the index of the end of the 
        complex. Indexes refer to the indices of x and yDat.
    """

    #Initialize empty list of bounds
    cBounds=[]

    #Iterate through cells of flux
    i=0
    while (i<len(x)):

        #Start tracking at a region that surpasses flux of edge
        if yDat[i]<complexLim:

            #Iterate through until reach next edge
            j=0
            while yDat[i+j]<complexLim: j=j+1

            #Check if the complex is big enough
            if j >minLength:

                #Check if there is enough absorption for the complex to
                #   be included
                cPeak = yDat[i:i+j].argmin()
                if yDat[cPeak+i]<fitLim:
                    cBounds.append([cPeak+i,i,i+j])

            i=i+j
        i=i+1

    i=0
    #Iterate through the bounds
    while i < len(cBounds):
        b=cBounds[i]

        #Check if the region needs to be divided
        if b[2]-b[1]>maxLength:

            #Find the minimum absorption in the middle two quartiles of
            #   the large complex
            q=(b[2]-b[1])/4
            cut = yDat[b[1]+q:b[2]-q].argmax()+b[1]+q

            #Only break it up if the minimum absorption is actually low enough
            if yDat[cut]>splitLim:

                #Get the new two peaks
                b1Peak = yDat[b[1]:cut].argmin()+b[1]
                b2Peak = yDat[cut:b[2]].argmin()+cut

                #add the two regions separately
                cBounds.insert(i+1,[b1Peak,b[1],cut])
                cBounds.insert(i+2,[b2Peak,cut,b[2]])

                #Remove the original region
                cBounds.pop(i)
                i=i+1
        i=i+1

    return cBounds

def _gen_flux_lines(x, linesP, speciesDict):
    """
    Calculates the normalized flux for a region of wavelength space
    generated by a set of absorption lines.

    Parameters
    ----------
    x : (N) ndarray
        Array of wavelength
    linesP: (3,) ndarray
        Array giving sets of line parameters in 
        form [[N1, b1, z1], ...]
    speciesDict : dictionary
        Dictionary containing all relevant parameters needed
        to create an absorption line of a given species (f,Gamma,lambda0)

    Returns
    -------
    flux : (N) ndarray
        Array of normalized flux generated by the line parameters
        given in linesP over the wavelength space given in x. Same size as x.
    """
    y=0
    for p in linesP:
        for i in range(speciesDict['numLines']):
            f=speciesDict['f'][i]
            g=speciesDict['Gamma'][i]
            wl=speciesDict['wavelength'][i]
            y = y+ _gen_tau(x,p,f,g,wl)
    flux = na.exp(-y)
    return flux

def _gen_tau(t, p, f, Gamma, lambda_unshifted):
    """This calculates a flux distribution for given parameters using the yt
    voigt profile generator"""
    N,b,z= p
    
    #Calculating quantities
    tau_o = 1.4973614E-15*N*f*lambda_unshifted/b
    a=7.95774715459E-15*Gamma*lambda_unshifted/b
    x=299792.458/b*(lambda_unshifted*(1+z)/t-1)
    
    H = na.zeros(len(x))
    H = voigt(a,x)
    
    tau = tau_o*H

    return tau

def _voigt_error(pTotal, x, yDat, yFit, speciesDict):
    """
    Gives the error of each point  used to optimize the fit of a group
        of absorption lines to a given flux profile.

        If the parameters are not in the acceptable range as defined
        in speciesDict, the first value of the error array will
        contain a large value (999), to prevent the optimizer from running
        into negative number problems.

    Parameters
    ----------
    pTotal : (3,) ndarray 
        Array with form [[N1, b1, z1], ...] 
    x : (N) ndarray
        array of wavelengths [nm]
    yDat : (N) ndarray
        desired normalized flux from fits of lines in wavelength
        space given by x
    yFit : (N) ndarray
        previous fit over the wavelength space given by x.
    speciesDict : dictionary
        dictionary containing all relevant parameters needed
        to create an absorption line of a given species (f,Gamma,lambda0)
        as well as max and min values for parameters to be fit

    Returns
    -------
    error : (N) ndarray
        the difference between the fit generated by the parameters
        given in pTotal multiplied by the previous fit and the desired
        flux profile, w/ first index modified appropriately for bad 
        parameter choices
    """

    pTotal.shape = (-1,3)
    yNewFit = _gen_flux_lines(x,pTotal,speciesDict)

    error = yDat-yFit*yNewFit
    error[0] = _check_params(pTotal,speciesDict)

    return error

def _check_params(p, speciesDict):
    """
    Check to see if any of the parameters in p fall outside the range 
        given in speciesDict.

    Parameters
    ----------
    p : (3,) ndarray
        array with form [[N1, b1, z1], ...] 
    speciesDict : dictionary
        dictionary with properties giving the max and min
        values appropriate for each parameter N,b, and z.

    Returns
    -------
    check : int
        0 if all values are fine
        999 if any values fall outside acceptable range
    """
    check = 0
    if any(p[:,0] > speciesDict['maxN']) or\
          any(p[:,0] < speciesDict['minN']) or\
          any(p[:,1] > speciesDict['maxb']) or\
          any(p[:,1] < speciesDict['minb']) or\
          any(p[:,2] > speciesDict['maxz']) or\
          any(p[:,2] < speciesDict['minz']):
              check = 999
    return check


def _output_fit(lineDic, file_name = 'spectrum_fit.h5'):
    """
    This function is designed to output the parameters of the series
    of lines used to fit an absorption spectrum. 

    The dataset contains entries in the form species/N, species/b
    species/z, and species/complex. The ith entry in each of the datasets
    is the fitted parameter for the ith line fitted to the spectrum for
    the given species. The species names come from the fitted line
    dictionary.

    Parameters
    ----------
    lineDic : dictionary
        Dictionary of dictionaries representing the fit lines. 
        Top level keys are the species given in orderFits and the corresponding
        entries are dictionaries with the keys 'N','b','z', and 'group#'. 
        Each of these corresponds to a list of the parameters for every
        accepted fitted line. 
    fileName : string, optional
        Name of the file to output fit to. Default = 'spectrum_fit.h5'

    """
    f = h5py.File(file_name, 'w')
    for ion, params in lineDic.iteritems():
        f.create_dataset("{0}/N".format(ion),data=params['N'])
        f.create_dataset("{0}/b".format(ion),data=params['b'])
        f.create_dataset("{0}/z".format(ion),data=params['z'])
        f.create_dataset("{0}/complex".format(ion),data=params['group#'])
    print 'Writing spectrum fit to {0}'.format(file_name)

