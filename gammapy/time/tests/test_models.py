# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import pytest
from numpy.testing import assert_allclose
from ..models_3 import PhaseCurve

time = 46601.0
reference = 43366.275 
phase0 = 0.0
f0 = 0.5
f1 = 0.0
f2 = 0.0

filename='/photon1/users/labsaha/gammapy_analysis/data/1dc/models/phasecurve_LSI_DC.fits'

def test_model():
    
    a = PhaseCurve(time,reference,phase0,f0,f1,f2)
    phase = a.phase(time)
    norm_const = a.evaluate(filename,time)
    assert_allclose(phase,0.36249 , rtol=1e-2)
    assert_allclose(norm_const,0.05 , rtol=1e-2)
