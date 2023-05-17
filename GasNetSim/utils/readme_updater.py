#   #!/usr/bin/env python
#   -*- coding: utf-8 -*-
#   ******************************************************************************
#     Copyright (c) 2023.
#     Developed by Yifei Lu
#     Last change on 5/16/23, 1:49 PM
#     Last change by yifei
#    *****************************************************************************

import re
from pathlib import Path

current_file_path = Path(__file__)

with open(current_file_path.parent.parent.parent / 'requirements.txt', 'r') as requirements_file:
    requirements = requirements_file.read()

version_cmp = ['<=', '<', '!=', '==', '>=', '>', '~=', '===']

for l in requirements.splitlines():
    package = re.findall(r'\b[a-zA-Z]+\b', l)
    if len(package) == 1:
        requirements = requirements.replace(package[0], f"- ``{package[0]}``")
    else:
        print(package)

with open(current_file_path.parent.parent.parent / 'README.md', 'r+') as readme_file:
    readme_contents = readme_file.read()

    # Find the section in the README.md file where you want to update the dependencies
    dependencies_section = re.search(r'<!-- Dependencies -->(.*?)<!-- End Dependencies -->', readme_contents, re.DOTALL)

    if dependencies_section:
        updated_dependencies_section = f"<!-- Dependencies -->\n{requirements}\n<!-- End Dependencies -->"
        updated_readme_contents = re.sub(r'<!-- Dependencies -->(.*?)<!-- End Dependencies -->',
                                         updated_dependencies_section, readme_contents, flags=re.DOTALL)

        # Move the file cursor to the beginning of the file and overwrite the contents
        readme_file.seek(0)
        readme_file.write(updated_readme_contents)
        readme_file.truncate()
