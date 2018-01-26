"""A setuptools based setup module for bemio
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

setup(

    name='python_scripts',

    version='beta',

    description='Michael Lawson\'s Python scripts',

    long_description='',

    url='https://github.com/lawsonro3/python_scripts',

    author='Michael Lawson',

    author_email='michael.lawson@nrel.gov',

    license='Apache 2.0',

    classifiers=[
        'Development Status :: Release',
        'Intended Audience :: Ocean research community',
        'License ::Apache 2.0',
        'Programming Language :: Python :: 3'
    ],

    keywords='python_scripts, Michael Lawson',

    packages=find_packages(exclude=['doc', 'tutorials']),

    install_requires=['numpy', 'scipy', 'matplotlib'],

    extras_require={
        'dev': [],
        'test': [],
    },

)
