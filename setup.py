import setuptools
import versioneer

with open("README.md", "r") as fh:
    long_description = fh.read()

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
    include_package_data=True,
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
