.. PEGG documentation master file, created by
   sphinx-quickstart on Wed Sep 28 12:28:06 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

|PEGG| PEGG: Prime Editing Guide Generator 
============================================

.. |PEGG| image:: PEGG_3.png
   :width: 200px
   :height: 200px

`Click here to read the bioRxiv preprint <https://www.biorxiv.org/content/10.1101/2022.10.26.513842v4>`_ 
******************************************************************************************************************

PEGG is a python module that designs prime editing guide RNAs (pegRNAs) for use in precision genome editing.
Unlike the existing, web-based programs for pegRNA design, PEGG is suitable for designing thousands of pegRNAs at once, giving users the ability to design entire libraries of pegRNAs.

PEGG's main functions are:

(1) Generating pegRNAs based on a list of input mutations

(2) Ranking and filtering these pegRNAs based on their properties

(3) Automated oligo generation (with the option to include a synthetic "sensor" region)

(4) Automated pegRNA library design with included silent substitution control pegRNAs

(5) Visualization Tools for individual pegRNA design and library design

Upcoming versions of PEGG will also include the ability to map mutations from the human to the mouse genome.

Installation
**************
PEGG is available through the python package index. To install, use pip: 

.. code-block:: python

   pip install pegg

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart
   jupyter_tutorial
   PEGG
   about

PEGG is an open source python package. If you use PEGG, please cite it using the following citation:

Gould, S.I., Sánchez-Rivera, F.J. (2022). PEGG: A computational pipeline for rapid design of prime editing guide RNAs and sensor libraries. *bioRxiv*. DOI: https://doi.org/10.1101/2022.10.26.513842.

Function Index
***************
* :ref:`genindex`
