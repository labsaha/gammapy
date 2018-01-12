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
 This is to calculate phase of a periodic system

# TODO: 1. Checking table values if the phases are in increasing order
       2. Check normalization constant to set maximum values
       
 parameters
 ----------
 time: event time or obervation time in MJD

 parameters : some parameters of the model need to be as input parameter

"""

class PhaseCurve(object):
      """ This is class docstring """
      
      def __init__(self,time,table,reference = 43366.275 ,phase0 = 0.0,f0 = 4.367575E-7,f1 = 0.0,f2 = 0.0):
          self.time = time
          self.table = table
          self.parameters = ParameterList([ 
          Parameter('reference', reference),
          Parameter('phase0', phase0),
          Parameter('f0', f0),
          Parameter('f1', f1),
          Parameter('f2', f2)])
     

      def phase(self,time):
          """ To assign the phase to the time """ 
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
     
      @classmethod
      def read(cls,filename, hdu=None):
          """Read `phase and normalization constant` table from file.

          Parameters
          ----------
          filename : str
            File name
          hdu : int or str
            HDU number or name

          Returns
          -------
          pointing_info : `PointingInfo`
            Pointing info object
          """
           
          filename = make_path(filename)

#          if hdu is None:
#             hdu = 'POINTING'
          table = Table.read(str(filename))
          return cls(table=table) 
#          return table
      def evaluate(self,time):

           phase = self.phase(time)
#           filename='/photon1/users/labsaha/gammapy_analysis/data/1dc/models/phasecurve_LSI_DC.fits'          
#           tbdata = len(self.table)
           ss = 'Table length: {}\n'.format(len(self.table))
           print(ss)
#           print('Tbdata',tbdata)
#           phase_col = tbdata['PHASE']
#           norm_col  = tbdata['NORM']
#            
#           x_values = phase_col
#           f_values = norm_col  

#           f_norm = np.interp(x = phase,xp= x_values,fp = f_values, period = 1)
           f_norm  = 1.0
           return f_norm

