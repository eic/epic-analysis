:orphan:

The documentation is automatically built after commit to master branch.

By default the documentation should be written in reStructuredText. Markdown files are acceptible also but have limited 
formatting possibilities compared to reStructuredText. `reStructuredText documentation <https://docutils.sourceforge.io/docs/user/rst/quickref.html>`_


Local Build
-----------


Requirments:
~~~~~~~~~~~~


- `Sphinx <http://www.sphinx-doc.org/en/master>`_ - Python documentation generator
- `Read the Docs Sphinx Theme <https://sphinx-rtd-theme.readthedocs.io/en/stable/>`_ - Theme for final output
- `recommonmark <https://github.com/miyakogi/m2r>`_ - Markdown to reStructuredText


.. code: bash

   pip install sphinx sphinx_rtd_theme recommonmark


Building
~~~~~~~~

.. code: bash

    pip install --upgrade sphinx-autobuild sphinx_rtd_theme recommonmark

    # from project root
    sphinx-autobuild docs docs/_build/html

    # from docs root
    sphinx-autobuild . _build/html



Read the docs
~~~~~~~~~~~~~

The documentation is available at https://g4e.readthedocs.io

It updates automatically when changes are pushed do the repo

