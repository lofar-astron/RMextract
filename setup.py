#from distutils.core import setup,Extension
from setuptools import setup, find_packages,Extension
setup(name='RMextract',
version='0.1',
ext_modules=[Extension('_EMM_Model', ['EMM/EMM_Model.cc','EMM/GeomagnetismLibrary.c','EMM/EMM_Model_wrap.cc'])],
extras_require = {'RMextract.LOFAR_TOOLS':['pyrap','lofar.parmdb']},
packages= ['RMextract','RMextract/LOFAR_TOOLS','EMM'],
package_data={'EMM':['*COF']},
)
