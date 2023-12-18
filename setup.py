from setuptools import setup
import setuptools

setup(
    include_package_data=True,
    name = 'pegg',
    author = 'Samuel Gould',
    author_email = 'samgould@mit.edu',
    url = 'https://github.com/samgould2/PEGG2.0',
    version = '2.0.9',
    description = 'Prime Editing Guide Generator',
    #package_dir = {'pegg': ''},
    #packages = setuptools.find_packages(), #['pegg'],
    #package_dir = {'pegg': '', 'azimuth': 'pegg/bin/Azimuth-2.0/azimuth'},
    packages = ['pegg'],
    #py_modules= ["pegg.prime", "pegg.base", "pegg.library", "pegg.crisporEffScores"],

    install_requires = ["Bio>=1.4.0",
        "cyvcf2>=0.30.18",
        "matplotlib>=3.5.1",
        "mock>=4.0.3",
        "numpy>=1.21.5",
        "pandas>=1.4.2",
        "seaborn>=0.11.2",
        "setuptools>=61.2.0",
        "Sphinx>=4.4.0",
        "scikit-learn==1.1.1",
        "regex>=2023.8.8",
        "importlib"
    ],

    classifiers=[
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.9",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent"
        ]
)