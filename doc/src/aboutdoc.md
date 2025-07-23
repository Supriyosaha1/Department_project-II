# About this documentation

The documentation is generated using git pages. 
The source files are in the directory ``doc/src/``. 

To contribute, you can edit these files, generate locally the doc, check it, and push. 

To build locally the documentation

```
pip install -r doc/sphinx_requirements.txt
```

```
sphinx-build -M html doc doc/_build
```

This will generate the documentation in the directory ``doc/_build/html``



