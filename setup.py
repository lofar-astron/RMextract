import os
from setuptools import find_packages
import numpy
from numpy.distutils.core import setup, Extension
from numpy.distutils.command.build_src import build_src as np_build_src
from numpy.distutils.misc_util import has_f_sources


def read(rel_path):
    """Function read() was copied from setup.py in Pip package."""
    here = os.path.abspath(os.path.dirname(__file__))
    # intentionally *not* adding an encoding option to open, See:
    #   https://github.com/pypa/virtualenv/issues/201#issuecomment-3145690
    with open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()


class build_src(np_build_src):
    """
    Override `build_src` class in `numpy.distutils.command.build_src`. This is necessary,
    because, when using `pip` (i.e. `setuptools`) to install from a source distribution,
    the files `fortranobject.h` and `fortranobject.` are not copied to the build source tree.
    However, these files are needed to build extension modules from Fortran sources. So we need
    to take care of that here, by overriding the `run` method in `build_src`.
    """

    def run(self):
        # Assume f2py source directory is `numpy/f2py/src`.
        srcdir = os.path.join(os.path.dirname(numpy.__file__), "f2py", "src")
        for ext in self.extensions:
            # If the extension contains Fortran sources, we need to copy `fortranobject.h` and
            # `fortranobject.c` from the f2py source directory to the extension's build directory.
            if has_f_sources(ext.sources):
                # The build source tree must have the same structure as the source tree. So the
                # name of the destination directory is the concatenation of the name of the build
                # source tree and the name of the directory containing the Fortran source files.
                # Here we assume that all the Fortran source files are in the same directory,
                # so we simply take the first one in the list to determine the source directory.
                destdir = os.path.join(self.build_src, os.path.dirname(ext.sources[0]))
                if not os.path.exists(destdir):
                    os.makedirs(destdir)
                for fn in ("fortranobject.h", "fortranobject.c"):
                    self.copy_file(os.path.join(srcdir, fn), os.path.join(destdir, fn))
        # Call the `run()` method of the parent
        np_build_src.run(self)


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
    version="0.4.1",
    url="https://github.com/lofar-astron/RMextract",
    project_urls={"Source": "https://github.com/lofar-astron/RMextract"},
    author="Maaike Mevius",
    author_email="mevius@astron.nl",
    description="Extract TEC, vTEC, Earthmagnetic field and Rotation Measures from GPS "
    "and WMM data for radio interferometry observations",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    cmdclass={"build_src": build_src},
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
