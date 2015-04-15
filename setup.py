from setuptools import setup, find_packages,Extension
import sys

packages=['RMextract','EMM']
scripts = ['RMextract/URL_download.py']
if "--add-lofar-utils" in sys.argv:
    packages.append("RMextract/LOFAR_TOOLS")
    scripts.append("RMextract/LOFAR_TOOLS/createRMParmdb")
    sys.argv.remove("--add-lofar-utils")


setup(name='RMextract',
version='0.1',
ext_modules=[Extension('_EMM_Model', ['EMM/EMM_Model.cc','EMM/GeomagnetismLibrary.c','EMM/EMM_Model_wrap.cc'])],
packages= packages,
package_data={'EMM':['*COF']},
 scripts = scripts
)
