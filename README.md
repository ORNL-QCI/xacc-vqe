[![Build Status](http://ci.eclipse.org/xacc/buildStatus/icon?job=xacc-vqe-ci)](http://ci.eclipse.org/xacc/job/xacc-vqe-ci/)

# XACC Plugins and Programs for the Variational Quantum Eigensolver

Installation
------------
With the XACC framework installed, users can choose a couple ways to install these plugins - using Python/Pip
```bash
$ python -m pip install --user .
```
or CMake and Make without Python support
```bash
$ mkdir build && cd build
$ cmake .. -DXACC_DIR=$HOME/.xacc (or wherever you installed XACC)
$ make install 
```
or with Python support
```bash
$ export PY_INC_DIR=$(python -c "import sysconfig; print(sysconfig.get_paths()['platinclude'])")
$ cmake .. -DXACC_DIR=$(python -m pyxacc -L) -DPYTHON_INCLUDE_DIR=$PY_INC_DIR
$ make install
```

Documentation
-------------

* [XACC Website and Documentation ](https://xacc.readthedocs.io)

Questions, Bug Reporting, and Issue Tracking
--------------------------------------------

Questions, bug reporting and issue tracking are provided by GitHub. Please
report all bugs by creating a new issue with the bug tag. You can ask
questions by creating a new issue with the question tag.

License
-------

XACC-VQE is licensed - [BSD 3-Clause](LICENSE).
