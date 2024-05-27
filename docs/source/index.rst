.. bc-pre-processing documentation master file, created by
   sphinx-quickstart on Thu May 23 17:19:03 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _pTatin3d: https://github.com/laetitialp/ptatin-gene
.. _PETSc: https://petsc.org

Welcome to bc-pre-processing's documentation!
=============================================
`bcpy` is a python module designed to evaluate symbolic (mathematic) 
expression to build analytical velocity functions varying in space and time, 
define the rheological parameters of a long-term geodynamic model, 
and generate input file for `pTatin3d`_.
It is a pre-processing tool.

This module can:

- evaluate and print mathematical expression for the velocity and initial plastic strain distribution 

- attribute rheological parameters to regions identified by a tag (integer value)

- handle simple mesh refinement for a structured mesh using linear interpolation

- generate options file for pTatin3d simulations

Some examples can be found in *scripts*  subdirectory.
Check out the :doc:`usage` section for further information, including how to
:ref:`install <installation>` the project.

.. note::
    This module is still under development.

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Contents
--------

.. toctree:: 

   usage
   bcpy_docs