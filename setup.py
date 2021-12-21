# TODO: Write a setup script for this project's packages

# Example:
"""
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="{{PROJECT-NAME}}",
    version="0.0.1",
    author="IEK-10, Forschungszentrum Jülich",
    author_email="author@fz-juelich.de",
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://jugit.fz-juelich.de/iek-10/sampleproject",
    project_urls={
        "Bug Tracker": "https://jugit.fz-juelich.de/iek-10/sampleproject/-/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "{{PROJECT-NAME}}"},
    packages=setuptools.find_packages(where="{{PROJECT-NAME}}"),
    python_requires=">=3.6",
)
"""
