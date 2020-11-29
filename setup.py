import setuptools
from setuptools import setup
setup(name='hzdplugins',
version='0.0.10',
description='plugins for my own research',
author='Zheng-Da He',
author_email='z.he@fz-juelich.de, jameshzd@mail.ustc.edu.cn, zhengdahe.electrocatalysis@gmail.com',
license='MIT',
packages=setuptools.find_packages(),
homepage='https://hzdplugins.readthedocs.io/en/latest/',
install_requires=['ase', 'pymatgen', 'aiida-core', 'aiida-quantumespresso'],
zip_safe=False)
