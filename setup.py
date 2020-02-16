import sys
from os.path import join
import setuptools
import versioneer
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy

COMMON_INCLUDE_DIRS = [
    './geoist/magmod',
    './geoist/magmod/include',
    join(sys.prefix, 'include'),
]
COMMON_INCLUDE_DIRS.append(numpy.get_include())

with open("README.md", "r") as fh:
    long_description = fh.read()

extensions = [
    Extension("geoist.pfm._prism",
              ["geoist/pfm/_prism.pyx"],
              include_dirs=[numpy.get_include()]
    ),
    Extension(
            'geoist.magmod._pymm',
            sources=[
                'geoist/magmod/pymm.c',
            ],
            libraries=[],
            library_dirs=[],
            include_dirs=COMMON_INCLUDE_DIRS,
        ),
    Extension(
            'geoist.magmod._pysunpos',
            sources=[
                'geoist/magmod/pysunpos.c',
            ],
            libraries=[],
            library_dirs=[],
            include_dirs=COMMON_INCLUDE_DIRS,
        ),
    Extension(
            'geoist.magmod._pytimeconv',
            sources=[
                'geoist/magmod/pytimeconv.c',
            ],
            libraries=[],
            library_dirs=[],
            include_dirs=COMMON_INCLUDE_DIRS,
        ),
]
install_requires = ["numpy","matplotlib","scipy","h5py","numba","pandas",
                    "pytest","future","Cython",'statsmodels>=0.9.0',"PyWavelets",
                    "seaborn","patsy", "appdirs"]
setuptools.setup(
    name="geoist",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="IGP-GRAVITY",
    author_email="gravity_igp@sina.com",
    description="An Open-Source Geophysical Python Library for Geoscience Prototype Research",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/igp-gravity/geoist",
    packages=setuptools.find_packages(),
    ext_modules=cythonize(extensions),
    include_package_data=True,
    install_requires=install_requires,
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
