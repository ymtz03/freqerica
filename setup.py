# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='freqerica',
    version='0.1.0',
    description='Hamiltonian dynamics simulator with quantum computer',
    long_description=readme,
    author='Yuta Matsuzawa',
    author_email='matsuzawa@theoc.kuchem.kyoto-u.ac.jp',
    url='https://github.com/ymtz03/freqerica',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')),
    package_data={'freqerica':['output/template/*']},
    install_requires=['numpy','scipy','sympy','openfermion','qulacs','pyscf'],
    test_suite='tests',
)
