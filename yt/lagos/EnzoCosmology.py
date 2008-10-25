import numpy as na

kmPerMpc = 3.08567758e19

class EnzoCosmology(object):
    def __init__(self, HubbleConstantNow = 71.0,
                 OmegaMatterNow = 0.27,
                 OmegaLambdaNow = 0.73,
                 OmegaCurvatureNow = 0.0,
                 InitialRedshift = 99.0):
        self.HubbleConstantNow = HubbleConstantNow
        self.OmegaMatterNow = OmegaMatterNow
        self.OmegaLambdaNow = OmegaLambdaNow
        self.OmegaCurvatureNow = OmegaCurvatureNow
        self.InitialRedshift = InitialRedshift
        self.InitialTime = self.ComputeTimeFromRedshift(self.InitialRedshift)
        self.TimeUnits = self.ComputeTimeUnits()

    def ComputeTimeUnits(self):
        "Taken from CosmologyGetUnits.C in Enzo."
        # Changed 2.52e17 to 2.52e19 because H_0 is in km/s/Mpc, 
        # instead of 100 km/s/Mpc.
        return 2.52e19 / na.sqrt(self.OmegaMatterNow) / \
            self.HubbleConstantNow / na.power(1 + self.InitialRedshift,1.5)

    def ComputeRedshiftFromTime(self,time):
        "Compute the redshift from time after the big bang.  This is based on Enzo's CosmologyComputeExpansionFactor.C, but altered to use physical units."

        OmegaCurvatureNow = 1.0 - self.OmegaMatterNow - self.OmegaLambdaNow

        OMEGA_TOLERANCE = 1e-5
        ETA_TOLERANCE = 1.0e-10

        # Convert the time to Time * H0.
 
        TimeHubble0 = time * self.HubbleConstantNow / kmPerMpc
 
        # 1) For a flat universe with OmegaMatterNow = 1, it's easy.
 
        if ((na.fabs(self.OmegaMatterNow-1) < OMEGA_TOLERANCE) and
            (self.OmegaLambdaNow < OMEGA_TOLERANCE)):
            a = na.power(time/self.InitialTime,2.0/3.0)
 
        # 2) For OmegaMatterNow < 1 and OmegaLambdaNow == 0 see
        #    Peebles 1993, eq. 13-3, 13-10.
        #    Actually, this is a little tricky since we must solve an equation
        #    of the form eta - na.sinh(eta) + x = 0..
 
        if ((self.OmegaMatterNow < 1) and 
            (self.OmegaLambdaNow < OMEGA_TOLERANCE)):
            x = 2*TimeHubble0*na.power(1.0 - self.OmegaMatterNow, 1.5) / \
                self.OmegaMatterNow;
 
            # Compute eta in a three step process, first from a third-order
            # Taylor expansion of the formula above, then use that in a fifth-order
            # approximation.  Then finally, iterate on the formula itself, solving for
            # eta.  This works well because parts 1 & 2 are an excellent approximation
            # when x is small and part 3 converges quickly when x is large. 
 
            eta = na.power(6*x,1.0/3.0)                # part 1
            eta = na.power(120*x/(20+eta*eta),1.0/3.0) # part 2
            for i in range(40):                      # part 3
                eta_old = eta
                eta = na.arcsinh(eta + x)
                if (na.fabs(eta-eta_old) < ETA_TOLERANCE): 
                    break
                if (i == 39):
                    print "No convergence after %d iterations." % i
 
            # Now use eta to compute the expansion factor (eq. 13-10, part 2).
 
            a = self.OmegaMatterNow/(2.0*(1.0 - self.OmegaMatterNow))*\
                (na.cosh(eta) - 1.0)

        # 3) For OmegaMatterNow > 1 and OmegaLambdaNow == 0, use sin/cos.
        #    Easy, but skip it for now.
 
        if ((self.OmegaMatterNow > 1) and 
            (self.OmegaLambdaNow < OMEGA_TOLERANCE)):
            print "Never implemented in Enzo, not implemented here."
            return 0
 
        # 4) For flat universe, with non-zero OmegaLambdaNow, see eq. 13-20.
 
        if ((na.fabs(OmegaCurvatureNow) < OMEGA_TOLERANCE) and
            (self.OmegaLambdaNow > OMEGA_TOLERANCE)):
            a = na.power(self.OmegaMatterNow / (1 - self.OmegaMatterNow),1.0/3.0) * \
                na.power(na.sinh(1.5 * na.sqrt(1.0 - self.OmegaMatterNow)*\
                                     TimeHubble0),2.0/3.0)


        redshift = (1.0/a) - 1.0

        return redshift

    def ComputeTimeFromRedshift(self,z):
        "Compute the time from redshift.  This is based on Enzo's CosmologyComputeTimeFromRedshift.C, but altered to use physical units."
        OmegaCurvatureNow = 1.0 - self.OmegaMatterNow - self.OmegaLambdaNow
 
        # 1) For a flat universe with OmegaMatterNow = 1, things are easy.
 
        if ((self.OmegaMatterNow == 1.0) and (self.OmegaLambdaNow == 0.0)):
            TimeHubble0 = 2.0/3.0/na.power(1+z,1.5)
 
        # 2) For OmegaMatterNow < 1 and OmegaLambdaNow == 0 see
        #    Peebles 1993, eq. 13-3, 13-10.
 
        if ((self.OmegaMatterNow < 1) and (self.OmegaLambdaNow == 0)):
            eta = na.arccosh(1 + 2*(1-self.OmegaMatterNow)/self.OmegaMatterNow/(1+z))
            TimeHubble0 = self.OmegaMatterNow/(2*na.power(1.0-self.OmegaMatterNow, 1.5))*\
                (na.sinh(eta) - eta)
 
        # 3) For OmegaMatterNow > 1 and OmegaLambdaNow == 0, use sin/cos.
 
        if ((self.OmegaMatterNow > 1) and (self.OmegaLambdaNow == 0)):
            eta = na.acos(1 - 2*(1-self.OmegaMatterNow)/self.OmegaMatterNow/(1+z))
            TimeHubble0 = self.OmegaMatterNow/(2*na.power(1.0-self.OmegaMatterNow, 1.5))*\
                (eta - na.sin(eta))
 
        # 4) For flat universe, with non-zero OmegaLambdaNow, see eq. 13-20.
 
        if ((na.fabs(OmegaCurvatureNow) < 1.0e-3) and (self.OmegaLambdaNow != 0)):
            TimeHubble0 = 2.0/3.0/na.sqrt(1-self.OmegaMatterNow)*\
                na.arcsinh(na.sqrt((1-self.OmegaMatterNow)/self.OmegaMatterNow)/ \
                               na.power(1+z,1.5))
  
        # Now convert from Time * H0 to time.
  
        time = TimeHubble0 / (self.HubbleConstantNow/kmPerMpc)
    
        return time
