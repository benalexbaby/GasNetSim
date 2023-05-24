![](docs/GasNetSim_Logo.svg)

[//]: # ([![PyPI]&#40;https://badge.fury.io/py/GasNetSim.svg&#41;]&#40;https://badge.fury.io/py/GasNetSim&#41;)
[![Python 3.9](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/release/python-370/)
[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
[![DOI:10.1109/OSMSES54027.2022.9769148](https://zenodo.org/badge/DOI/10.1109/OSMSES54027.2022.9769148.svg)](https://doi.org/10.1109/OSMSES54027.2022.9769148)

# **GasNetSim**

*GasNetSim* is a simulation package designed for gas network steady-state simulation. 
It supports the steady-state natural gas network simulations with different gas mixture
compositions, thus enabling accurate analysis of the impacts of hydrogen injection on the gas network.  
Moreover, users have the flexibility to modify this tool and implement their own desired
gas mixture modeling approaches.
Future work will be carried out to include gas storage units and to take into account 
the dynamic behavior of the gas network so that short-term simulations can be performed.

## Installation
Currently, it is only supported using source files.

## License

The project is released under the terms of the [MPL 2.0](https://mozilla.org/MPL/2.0/).

## Dependencies
<!-- Dependencies -->
- ``numpy``>=1.19.2
- ``matplotlib``>=3.3.2
- ``scipy``>=1.5.2
- ``pandas``>=1.1.3
- ``pytest``>=6.2.5
- ``fluids``>=0.1.86
- ``pint``>=0.18
- ``setuptools``>=60.9.3
- ``requests``>=2.25.1
- ``pyparsing``~=3.0.7
- ``cantera``~=2.6.0
- ``thermo``~=0.2.23
- ``tqdm``>=4.64.1
- ``seaborn``>=0.12.2
- ``networkx``~=3.1
<!-- End Dependencies -->

For the ``thermo`` package, the version used in this repo is 0.1.40. Because there are some changes 
and new features included in the newer versions. The source files of the `thermo` package is directly
stored in this repo. It will be updated in the future.

[//]: # (## Discussion)

[//]: # ()
[//]: # (You can connect with the community in a variety of ways...)

[//]: # ()
[//]: # (- [Mailing list]&#40;https://lists.lfenergy.org/g/xxxx-discussion&#41;)

[//]: # (- [#{{**PROJECT-NAME**}} channel on LF Energy Slack]&#40;https://slack.lfenergy.org&#41;)

[//]: # (- Other communication channels, e.g. Discord, Slack, Skype, Mattermost, FZJ Rocket Chat, ...)

[//]: # (## Contributing)

[//]: # ()
[//]: # (_**TODO** Provide contributing guidelines here or point to a_)

[//]: # (_[CONTRIBUTING.md]&#40;CONTRIBUTING.md&#41; file if the contributing guidelines require_)

[//]: # (_more than just a few lines._)

## Reporting Issues

To report a problem, you can open an 
[issue](https://jugit.fz-juelich.de/iek-10/public/simulation/gasnetsim/-/issues)
in repository against a specific workflow. If the issue is sensitive in nature or
a security related issue, please do not report in the issue tracker but instead
email [Yifei Lu](yi.lu@fz-juelich.de).
