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

PEGG is a python package that designs prime editing guide RNAs (pegRNAs) and base editing guide RNAs (gRNAs) for use in precision genome editing.
Unlike the existing, web-based programs for pegRNA design, PEGG is suitable for designing thousands of pegRNAs at once, giving users the ability to design entire libraries of pegRNAs
and gRNAs. Uniquely, PEGG can design paired pegRNA or gRNA-sensor cassettes that include a synthetic version of the target locus, allowing for 
the calibration of guide editing activity in pooled screens (see above bioRxiv preprint for more information).

PEGG's main functions are:

(1) Generating pegRNAs or gRNAs based on a list of input mutations. The input format is extremely flexible, allowing for users to input a list of genome coordinates or sequences with desired edits.

(2) Ranking and filtering these pegRNAs based on their properties, including On-Target (Azimuth) Scores.

(3) Automated oligo generation (with the option to include a synthetic "sensor" region).

(4) Automated pegRNA/gRNA library design with included safe-targeting mutations, non-targeting guides, and silent substitution controls.

(5) Visualization Tools for pegRNA and gRNA design.

PEGG has recently been updated to version 2.0, with new features including (1) increased input mutation format flexibility,
(2) Dynamically computed Azimuth on-target scores, (3) a new base editing module, (4) improved library design functionality, as well as some minor bug fixes with INS/DEL design.

Installation
**************
PEGG is available through the python package index. To install, use pip: 

.. code-block:: python

   pip install pegg

Note
*****
PEGG has been tested with python version 3.9 and 3.10. To get it to install, you may need to use a `virtual environment <https://saturncloud.io/blog/how-to-install-python-39-with-conda-a-guide-for-data-scientists/>`_ :

.. code-block:: python

   conda create -n myenv python=3.9

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart
   jupyter_tutorial
   PEGG
   about
   other_design_tools

PEGG is an open source python package. If you use PEGG, please cite it using the following citation:

Gould, S. I., Wuest, A. N., Dong, K., Johnson, G. A., Hsu, A., Narendra, V. K., Levine, S. S., Liu, D. R., & Sánchez Rivera, F. J. (2022). High throughput evaluation of genetic variants with prime editing sensor libraries. bioRxiv. https://doi.org/10.1101/2022.10.26.513842

Function Index
***************
* :ref:`genindex`
