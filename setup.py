from setuptools import setup, find_packages

setup(
    name="sccompPy",  # Your package name
    version="0.1.0",  # Package version
    author="Chen Zhan",  # Your name or the author's name
    author_email="chen.zhan@adelaide.edu.au",  # Author's email
    description="A Python package for sccomp - Tests differences in cell type proportions and variability from single-cell data",  # Short description
    long_description=open("README.md").read(),  # Load README as the long description
    long_description_content_type="text/markdown",  # Set content type for long description
    url="https://github.com/MangiolaLaboratory/sccompPy",  # URL to the project (GitHub or another platform)
    packages=find_packages(),  # Automatically find and include all packages
    include_package_data=True,  # Include non-Python files specified in MANIFEST.in
    package_data={
        "sccompPy": ["data/*", "stan/*.stan"],  # Include all files in the data folder
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="GPLv3",  # Specify the license
    python_requires=">=3.7",  # Minimum Python version required
    install_requires=[
        "pandas>=1.0.0",        # For data manipulation
        "numpy>=1.18.0",        # For numerical computations
        "patsy>=0.5.0",         # For creating design matrices
        "cmdstanpy>=1.0.0",     # For Stan models
    ],
    extras_require={
        "dev": ["pytest", "black", "flake8"],  # Optional development dependencies
    },
)
