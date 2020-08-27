from setuptools import find_packages, setup

with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='espsim',
    version='0.0.1',
    author='Esther Heid',
    author_email='eheid@mit.edu',
    description='Scoring of shape and ESP similarity with RDKit',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/hesther/espsim',
    license='MIT',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent'
    ],
    keywords=[
        'chemistry',
        'electrostatic potential',
        'shape',
        'similarity',
        'RDKit'
    ]
)
