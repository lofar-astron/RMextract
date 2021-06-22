# import sys
from setuptools import find_packages
import numpy
from numpy.distutils.core import setup, Extension

packages = ['RMextract', 'RMextract.EMM']  # , 'EMM']
packages = find_packages(exclude=['RMextract/LOFAR_TOOLS'])

# scripts = ['RMextract/URL_download.py']
scripts = []
ext = [Extension('_EMM_Model',
                 ['RMextract/EMM/EMM_Model.cc', 'RMextract/EMM/GeomagnetismLibrary.c',
                  'RMextract/EMM/EMM_Model_wrap.cc'],
                 extra_compile_args=['-Wno-format-security'])
       ]

if False:  # "--add-lofar-utils" in sys.argv:
    packages.append("RMextract/LOFAR_TOOLS")
    scripts.append("RMextract/LOFAR_TOOLS/createRMParmdb")
    scripts.append("RMextract/LOFAR_TOOLS/createRMh5parm.py")
    scripts.append("RMextract/LOFAR_TOOLS/download_IONEX.py")
    # sys.argv.remove("--add-lofar-utils")

if True:  # "--add-iri" in sys.argv:
    packages.append("RMextract/pyiri")
    ext.append(Extension('_iri',
                         sources=['RMextract/pyiri/iri.pyf', 'RMextract/pyiri/cira.for', 'RMextract/pyiri/igrf.for', 'RMextract/pyiri/iridreg.for',
                                  'RMextract/pyiri/iriflip.for', 'RMextract/pyiri/irifun.for', 'RMextract/pyiri/irisub.for',
                                  'RMextract/pyiri/iritec.for', 'RMextract/pyiri/iriget.for'],
                         include_dirs=[numpy.get_include()])
               )
    packages.append("RMextract/pyiriplas")
    ext.append(Extension('_iriplas',
                         sources=['RMextract/pyiriplas/iriplas.pyf', 'RMextract/pyiriplas/igrf.for', 'RMextract/pyiriplas/irif2019.for',
                                  'RMextract/pyiriplas/iriplas_main.for', 'RMextract/pyiriplas/Iris2017.for',
                                  'RMextract/pyiriplas/indx2017.for'],
                         include_dirs=[numpy.get_include()])
               )
    # sys.argv.remove("--add-iri")

setup(name='RMextract',
      version='0.4',
      author='Maaike Mevius',
      author_email='mevius@astron.nl',
      url='https://github.com/lofar-astron/RMextract',
      ext_modules=ext,
      packages=packages,
      install_requires=['numpy', 'scipy', 'python-casacore'],
      package_data={'RMextract/EMM': ['*COF'],
                    'RMextract/pyiri': ['*dat', '*asc'],
                    'RMextract/pyiriplas': ['*dat', '*asc', 'kp*', '*ASC']},
      scripts=scripts
      )
