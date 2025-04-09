.. _m4opt:
---

1. Setting Up M4OPT
-------------------

Create and Activate a Conda Environment
---------------------------------------

```bash
conda create --name m4opt_env python=3.11
conda activate m4opt_env
```

Clone and Install M4OPT
-------------------
```bash
git clone https://github.com/m4opt/m4opt.git
cd m4opt/
pip install -e .
cd ..
```

Verify Installation
-------------------

```bash
m4opt --help
m4opt schedule --help
m4opt schedule --mission ultrasat
```

2. Install CPLEX
----------------

[Install CPLEX](https://m4opt.readthedocs.io/en/latest/install/cplex.html)


3. Installing Additional Dependencies

```bash
pip install joblib
pip install dask
```

---

4. Verifying Installation
-------------------------

```bash
m4opt schedule --help
```

For additional documentation, visit [M4OPT ReadTheDocs](https://m4opt.readthedocs.io/en/latest/install/index.html).
