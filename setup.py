from setuptools import find_packages,Extension
from numpy.distutils.core import setup
from numpy.distutils.core import Extension as Extension2
import sys

packages=['RMextract','EMM']
#scripts = ['RMextract/URL_download.py']
scripts = []
ext=[Extension2('_EMM_Model', ['EMM/EMM_Model.cc','EMM/GeomagnetismLibrary.c','EMM/EMM_Model_wrap.cc'],extra_compile_args=['-Wno-format-security'])]

if "--add-lofar-utils" in sys.argv:
    packages.append("RMextract/LOFAR_TOOLS")
    scripts.append("RMextract/LOFAR_TOOLS/createRMParmdb")
    scripts.append("RMextract/LOFAR_TOOLS/createRMh5parm.py")
    sys.argv.remove("--add-lofar-utils")
if "--add-iri" in sys.argv:
    packages.append("pyiri")
    ext.append(Extension2('_iri',sources=['pyiri/iri.pyf','pyiri/cira.for',  'pyiri/igrf.for',  'pyiri/iridreg.for',  'pyiri/iriflip.for',  'pyiri/irifun.for',  'pyiri/irisub.for',  'pyiri/iritec.for',  'pyiri/iriget.for']))
    sys.argv.remove("--add-iri")


setup(name='RMextract',
version='0.1',
      ext_modules=ext,
      packages= packages,
      package_data={'EMM':['*COF'],
                    'pyiri':['*dat','*asc']},
 scripts = scripts
)
