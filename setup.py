import setuptools
import versioneer
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy

with open("README.md", "r") as fh:
    long_description = fh.read()

extensions = [
    Extension("geoist.pfm._prism",
              ["geoist/pfm/_prism.pyx"],
              include_dirs=[numpy.get_include()]
    )
]
install_requires = ["numba","pandas","pytest","future","Cython"]
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
