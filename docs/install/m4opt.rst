.. _m4opt:

M4OPT Setup Guide
=================

1. Setting Up M4OPT
-------------------

Create and Activate a Conda Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   conda create --name m4opt_env python=3.11
   conda activate m4opt_env

Clone and Install M4OPT
~~~~~~~~~~~~~~~~~~~~~~~

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

`M4OPT Documentation <https://m4opt.readthedocs.io/en/latest/install/index.html>`_
