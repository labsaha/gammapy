from __future__ import absolute_import, division, print_function, unicode_literals
import logging
import numpy as np
import astropy
from astropy.table import Table
import astropy.units as u
from  ..utils.modeling import Parameter, ParameterList
from ..utils.scripts import make_path

__all__ = [
    'PhaseCurve',
]

log = logging.getLogger(__name__)

"""
 This is to calculate phase of a periodic system and to provide the nomralization factor for the phase using cirular interpolation. The required data for interpolation is read from a fits file

# TODO: 1. Checking table values if the phases are in increasing order
        2. Check normalization constant to set maximum values
       
 Parameters
 ----------
 time: event time or obervation time in MJD
 reference: The MJD value where phase is considered as 0.
 phase0: phase at the refernce MJD 
 f0:1st order derivative of the function phi with time 
 f1:2nd order derivative of the function phi with time  
 f2:3rd order derivative of the function phi with time 
 phase as function of time is computed using
 Φ(t)=Φ0+f0(t−t0)+(1/2)f1(t−t0)^2+(1/6)f2(t−t0)^3 
 filename = Name of the Fits file for 'PHASE' vs 'NORM' along with the path shold be given 
 
 Returns
 --------
 Phase corresponding to the time provided as input
 It also calculates the normalization constant for that phase

 Examples
 --------
 This shows how to get the phase and the normalization constant

    a = PhaseCurve(time,reference,phase0,f0,f1,f2)
    phase = a.phase(time)
    norm_const = a.evaluate(filename,time)

"""

class PhaseCurve(object):
      """ This is class docstring """
      
      def __init__(self,time,reference = 43366.275 ,phase0 = 0.0,f0 = 4.367575E-7,f1 = 0.0,f2 = 0.0):
          self.time = time
#          self.table = table
          self.parameters = ParameterList([ 
          Parameter('reference', reference),
          Parameter('phase0', phase0),
          Parameter('f0', f0),
          Parameter('f1', f1),
          Parameter('f2', f2)])
     

      def phase(self,time):
          """ To assign the phase to the time.
              It returns phase
          """
          c1 = 0.5
          c2 = 1.0/6.0
          pars = self.parameters
          reference = pars['reference'].value
          phase0    = pars['phase0'].value
          f0        = pars['f0'].value
          f1        = pars['f1'].value
          f2        = pars['f2'].value
           
          t = self.time - reference
          
          phase = phase0 + t * (f0 + t * (c1 * f1 +  c2 * f2 * t))
          
          phase = np.remainder(phase,1)
          return phase
     
      def read(self,filename, hdu=None):
          """Read `phase and normalization constant` table from file.

          Parameters
          ----------
          filename : str
            File name

          Returns
          -------
          an astropy table
          """
          table = Table.read(str(filename))
          return table

      def evaluate(self,filename,time):
           """
            Evaluates the normalization constant for the given phase using
            circular interpolation
           
           Parameters
           -----------
           filename: str
           time    : float
          
           returns
           -----------
           normalization constant
            
           """
           phase = self.phase(time)
           
           tbdata = self.read(filename)
           phase_col = tbdata['PHASE']
           norm_col  = tbdata['NORM']
           
           x_values = phase_col
           f_values = norm_col  
            
           f_norm = np.interp(x = phase,xp= x_values,fp = f_values, period = 1)

           return f_norm

