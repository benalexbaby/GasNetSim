# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Support for parallel pipelines
- Test for pipeline outlet temperature calculation function
- Test and validation with original [C++ implementation](https://github.com/usnistgov/AGA8/tree/master/AGA8CODE/C) of GERG-2008 EOS
- Added support of heating value calculation for GERG2008 gas mixture class
- Tests and validations for the implementation of the heating values calculation
- Tests and validations of gas mixture properties calculation using the _thermo_ package and the GERG-2008 EOS

### Changed
- Change `Cantera` dependency to `v3.0.0`
- Fixed pipeline outlet temperature calculation
- Fixed GERG-2008 implementation
- Set default gas mixture EOS for simulation to GERG-2008

### Removed
- Local adapted implementation of `thermo` package