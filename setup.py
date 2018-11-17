import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="geoist",
    version="0.0.1",
    author="IGP-GRAVITY",
    author_email="gravity_igp@sina.com",
    description="An Open-Source Geophysical Python Library for Geoscience Prototype Research",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/igp-gravity/geoist",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
