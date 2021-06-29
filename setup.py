import os
from setuptools import find_packages
import numpy
from numpy.distutils.core import setup, Extension

packages = find_packages(exclude=["RMextract.LOFAR_TOOLS"])

ext_modules = []
scripts = []

ext_modules.append(
    Extension(
        "RMextract.EMM._EMM_Model",
        sources=[
            os.path.join("RMextract", "EMM", f)
            for f in ("EMM_Model.cc", "GeomagnetismLibrary.c", "EMM_Model_wrap.cc")
        ],
        extra_compile_args=["-Wno-format-security"],
    )
)

ext_modules.append(
    Extension(
        "RMextract.pyiri._iri",
        sources=[
            os.path.join("RMextract", "pyiri", f)
            for f in (
                "iri.pyf",
                "cira.for",
                "igrf.for",
                "iridreg.for",
                "iriflip.for",
                "irifun.for",
                "irisub.for",
                "iritec.for",
                "iriget.for",
            )
        ],
        include_dirs=[numpy.get_include()],
    )
)

ext_modules.append(
    Extension(
        "RMextract.pyiriplas._iriplas",
        sources=[
            os.path.join("RMextract", "pyiriplas", f)
            for f in (
                "iriplas.pyf",
                "igrf.for",
                "irif2019.for",
                "iriplas_main.for",
                "Iris2017.for",
                "indx2017.for",
            )
        ],
        include_dirs=[numpy.get_include()],
    )
)

if "RMextract.LOFAR_TOOLS" in packages:
    scripts.extend(
        [
            os.path.join("RMextract", "LOFAR_TOOLS", f)
            for f in ("createRMParmdb", "createRMh5parm.py", "download_IONEX.py")
        ]
    )

setup(
    name="RMextract",
    version="0.4",
    author="Maaike Mevius",
    author_email="mevius@astron.nl",
    url="https://github.com/lofar-astron/RMextract",
    ext_modules=ext_modules,
    packages=packages,
    install_requires=["numpy", "scipy", "python-casacore"],
    package_data={
        "RMextract.EMM": ["*.COF"],
        "RMextract.pyiri": ["*.dat", "*.asc"],
        "RMextract.pyiriplas": ["*.dat", "*.asc", "kp*", "*.ASC"],
    },
    scripts=scripts,
)
