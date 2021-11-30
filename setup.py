import os
import os.path
from setuptools import find_packages, setup
import codecs

# this comes form here https://packaging.python.org/guides/single-sourcing-package-version/
def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

setup(
    #package_data={'': package_files('polyply/data')
    #              + package_files('polyply/tests/test_data'),},
    scripts=['bin/polyply_gui'],
    pbr=False,
    version=get_version("polyply_gui/__init__.py")
)
