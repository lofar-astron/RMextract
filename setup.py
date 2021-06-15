from numpy.distutils.core import setup
from numpy.distutils.core import Extension
import numpy
import sys

packages = ['RMextract', 'EMM']
# scripts = ['RMextract/URL_download.py']
scripts = []
ext = [Extension('_EMM_Model',
                 ['EMM/EMM_Model.cc', 'EMM/GeomagnetismLibrary.c', 'EMM/EMM_Model_wrap.cc'],
                 extra_compile_args=['-Wno-format-security'])
       ]

if "--add-lofar-utils" in sys.argv:
    packages.append("RMextract/LOFAR_TOOLS")
    scripts.append("RMextract/LOFAR_TOOLS/createRMParmdb")
    scripts.append("RMextract/LOFAR_TOOLS/createRMh5parm.py")
    scripts.append("RMextract/LOFAR_TOOLS/download_IONEX.py")
    sys.argv.remove("--add-lofar-utils")

if "--add-iri" in sys.argv:
    packages.append("pyiri")
    ext.append(Extension('_iri',
                         sources=['pyiri/iri.pyf', 'pyiri/cira.for', 'pyiri/igrf.for', 'pyiri/iridreg.for',
                                  'pyiri/iriflip.for', 'pyiri/irifun.for', 'pyiri/irisub.for',
                                  'pyiri/iritec.for', 'pyiri/iriget.for'],
                         include_dirs=[numpy.get_include()])
               )
    packages.append("pyiriplas")
    ext.append(Extension('_iriplas',
                         sources=['pyiriplas/iriplas.pyf', 'pyiriplas/igrf.for', 'pyiriplas/irif2019.for',
                                  'pyiriplas/iriplas_main.for', 'pyiriplas/Iris2017.for',
                                  'pyiriplas/indx2017.for'],
                         include_dirs=[numpy.get_include()])
               )
    sys.argv.remove("--add-iri")

setup(name='RMextract',
      version='0.4',
      ext_modules=ext,
      packages=packages,
      install_requires=['numpy', 'scipy', 'python-casacore'],
      package_data={'EMM': ['*COF'],
                    'pyiri': ['*dat', '*asc'],
                    'pyiriplas': ['*dat', '*asc', 'kp*', '*ASC']},
      scripts=scripts
      )
