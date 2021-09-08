from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.readlines()

long_description = """Functional profiling tool for 16S amplicon-based sequence
data. Based on SecondGenome's Piphillin webserver, since discontinued."""

setup(
    name = 'pyphillin',
    version = '1.0.0',
    author = 'Daniel McSkimming',
    author_email = 'dmcskimming@usf.edu',
    url = 'https://github.com/dmcskim/pyphillin',
    description = 'Functional profiling tool for 16S amplicon-based sequence\
    data.',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    license = 'GPLv3',
    packages = find_packages(),
    entry_points = {
        'console_scripts' : [ 'pyphillin = pyphillin.pyphillin:main' ]
    },
    classifiers = ("Programming Language :: Python :: 3",
                   "Operating System :: OS Independent",),
    install_requires = requirements,
    include_package_data = True,
    package_data = {
        '' : ['data/*']
    },
)
