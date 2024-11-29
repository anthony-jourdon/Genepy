from setuptools import setup, find_packages

setup(
    name='genepy',
    version='1.1.3',
    packages=find_packages(include=['genepy', 'genepy.*']),
    include_package_data=True,
    package_data={
        'genepy.material_params': ['arrhenius_flow_laws.json'],
    },
    install_requires=[
        # List your package dependencies here
        "numpy",
        "sympy",
        "pyvista",
        "matplotlib",
        "gmsh"
    ],
    entry_points={
        'console_scripts': [
            # Define command-line scripts here
        ],
    },
    author='Anthony Jourdon',
    author_email='',
    description='Python package to build input files for the geodynamic code pTatin3d',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/anthony-jourdon/Genepy',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.12',
)