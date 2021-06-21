The NiBabies Docker/Singularity wrapper
---------------------------------------

NiBabies is a functional magnetic resonance image pre-processing pipeline
optimized for infant and neonate MRI. It is designed to provide an easily
accessible, state-of-the-art interface that is robust to differences in
scan acquisition protocols and that requires minimal user input, while
providing easily interpretable and comprehensive error and output reporting.

This is a Python wrapper to run NiBabies.
It generates the appropriate Docker or Singularity commands, providing an
intuitive interface to running the fMRIPrep workflow in whichever environment.
Docker or Singularity must be installed, and in the case of Docker, running.
Installations can be check by running ::

  docker info  # Docker
  singularity version  # Singularity

Please report any feedback to our `GitHub repository
<https://github.com/nipreps/nibabies>`_ and do not
forget to `credit <https://fmriprep.readthedocs.io/en/latest/citing.html>`_ all
the authors of software that NiBabies uses.


Usage
-----

Example Docker usage ::

  nibabies-wrapper docker <data-path> <output-path> participant <nibabies-arguments>

Example Singularity usage ::

  nibabies-wrapper singularity <data-path> <output-path> participant -i <img-path> <nibabies-arguments>
