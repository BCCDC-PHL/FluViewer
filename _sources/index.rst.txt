.. BCCDC-PHL TB Genomics Database documentation master file, created by
   sphinx-quickstart on Thu Dec 15 16:15:13 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

BCCDC-PHL FluViewer
==========================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

fluviewer.fluviewer
===================
This module controls the overall flow of analysis by running each analysis stage.

It runs the following analysis stages:

0. Read Normalization
1. Assemble Contigs
2. BLAST Contigs
3. Scaffolding
4. Read Mapping
5. Variant Calling
6. Consensus Calling
7. Summary Reporting

.. automodule:: fluviewer.fluviewer
   :members:

fluviewer.analysis
==================
This module includes the steps required for each analysis stage.

.. automodule:: fluviewer.analysis
   :members:

fluviewer.database
==================
This module includes functions used to check the integrity of the
FluViewer database prior to analysis.
      
.. automodule:: fluviewer.database
   :members:

fluviewer.parsers
=================
This module includes functions related to parsing output files.
      
.. automodule:: fluviewer.parsers
   :members:

fluviewer.report
================
This module includes functions related to building the final summary report.


.. automodule:: fluviewer.report
   :members:

fluviewer.plots
===============
This module includes functions related to generating plots to visualize FluViewer outputs.


.. automodule:: fluviewer.plots
   :members:


fluviewer.cli_args
===============
This module includes functions related to parsing and validating command-line arguments.


.. automodule:: fluviewer.cli_args
   :members:


fluviewer.logging_config
===============
This module includes functions related to logging configuration.


.. automodule:: fluviewer.logging_config
   :members:

