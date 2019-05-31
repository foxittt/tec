
"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'DESCRIPTION.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='tec-embrace',
    version='0.1',
    url='https://github.com/embrace-inpe/tec',
    license='MIT',
    author='Estudo e Monitoramento Brasileiro do Clima Espacial (EMBRACE/INPE)',
    author_email='desenvolvimento.embrace@gmail.com',
    description='TEC and receiver bias estimation model',
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    long_description=open('README.md').read(),
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
    ],
    zip_safe=False)
