.. _m4opt:

.. ..:math:`\mathrm{M^4OPT}` Setup
.. =============

1. Setting Up |M4OPT|
---------------------

Clone and install |M4OPT|
~~~~~~~~~~~~~~~~~~~~~~~~~

The recommended way to install `M⁴OPT <https://m4opt.readthedocs.io/en/latest/install/index.html>`_ is using pip:

.. code-block:: bash

   pip install m4opt

Alternatively, you can clone the repository and install it in editable mode
(if you want to explore or modify the source code):

.. code-block:: bash

   git clone https://github.com/m4opt/m4opt.git
   cd m4opt/
   pip install -e .
   cd ..

Verify Installation
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   m4opt --help
   m4opt schedule --help
   m4opt schedule --mission ultrasat

2. Install CPLEX
----------------

To install CPLEX, follow the instructions at:

`Install CPLEX <https://m4opt.readthedocs.io/en/latest/install/cplex.html>`_

3. Installing Additional Dependencies
-------------------------------------

.. code-block:: bash

   pip install joblib
   pip install dask

4. Verifying Installation
-------------------------

.. code-block:: bash

   m4opt schedule --help

For additional documentation:

`M4OPTocumentation <https://m4opt.readthedocs.io/en/latest/install/index.html>`_
