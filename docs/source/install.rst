Installation
============

This section provide guidelines for installing XACC-VQE and its TPLs.

Install third-party libraries
-----------------------------

First, see `XACC Install
<http://xacc.readthedocs.io/en/latest/install.html>`_ to install XACC. 

The following third party libraries (TPLs) are used by XACC-VQE:

+------------------------+------------+-----------+
| Packages               | Dependency | Version   |
+========================+============+===========+
| MPI                    | Required   | See below |
+------------------------+------------+-----------+
| Others...              | Required   | 1.59.0+   |
+------------------------+------------+-----------+

Note that you must have a C++11 compliant compiler. 
For GCC, this means version 4.8.1+, for Clang, this means 3.3+.

These dependencies are relatively easy to install on popular operating
systems. Any of the following commands will work (showing with and without MPI):

.. code::

   $ (macosx) brew install boost
   $ (fedora) dnf install boost-devel
   $ (ubuntu) apt-get install libboost-all-dev

Build XACC-VQE
-----------

Clone the XACC-VQE repository:

.. code::

   $ git clone https://github.com/ornl-qci/xacc-vqe

XACC-VQE requires CMake 3.2+ to build. Run the following to
configure and build XACC-VQE:

.. code:: bash

   $ cd xacc-vqe && mkdir build && cd build
   $ cmake ..
   $ make install # can pass -jN for N = number of threads to use

This will build, test, and install XACC-VQE to ``/usr/local/xacc``
(Pass ``-DCMAKE_INSTALL_PREFIX=$YOURINSTALLPATH`` to install it somewhere else).

