# RASCAS installation

Before [downloading](#download-the-rascas-code) and [building](#building-the-code) the RASCAS code, make sure to have installed on your laptop the following apps/softwares. 

## Environment

To compile and run the RASCAS code, one needs to have:
- [ ] a Fortran compiler, for instance gfortran, the [GNU Fortran compiler](https://gcc.gnu.org/fortran/)
- [ ] the [GNU make tool](https://www.gnu.org/software/make/)
- [ ] an MPI library, for instance [OpenMPI](https://www.open-mpi.org/)

To visualise the results and to run the notebooks tutorials, one needs:
- [ ] a minimal Python installation with python3, numpy, scipy, matplotlib, and astropy
- [ ] the [Jupyter notebook](https://jupyter.org/)

It is also very useful to use git (but not mandatory). 

### expert mode

If you are familiar with the installation of packages, you can install these packages with your favorite package manager. 

### non-expert mode

Alternatively, follow the instructions below to set up a RASCAS environment using [conda](https://conda.io/projects/conda/en/latest/index.html). Conda is available for Windows, Linux, and macOS. 

1. Install miniconda following the [installation guide](https://docs.conda.io/en/latest/miniconda.html), or if you already have a version of conda installed, make sure to have an up-to-date version by running the following command in a terminal.

```
conda update -n base conda
```

2. Download the environment file [rascas-environment.yml](https://git-cral.univ-lyon1.fr/rascas/rascas/-/blob/master/doc/rascas-environment.yml)

3. Create the environment from the ``rascas-environment.yml`` file by executing the command below in a terminal.
```
conda env create -f rascas-environment.yml
```

4. Activate the new environment by executing the command below in a terminal.
```
conda activate rascas-env
```

5. Verify that the new environment was installed correctly with
```
conda env list
```

A detailed documentation about conda environments is available [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#)


## Download the RASCAS code

There are two options here. The first one (recommended) is to get the code using git. If you prefer not to use git, you can download a tarball of the code. 

### Using git
In a terminal, run the following commands. 
```
git clone https://git-cral.univ-lyon1.fr/rascas/rascas.git
cd rascas
```

### Simple download

Alternatively, you can obtain a tarball of the code using either `curl` or `wget`
```
curl -O https://git-cral.univ-lyon1.fr/rascas/rascas/-/archive/master/rascas-master.tar.gz
```
or
```
wget https://git-cral.univ-lyon1.fr/rascas/rascas/-/archive/master/rascas-master.tar.gz
```
Then, you can untar and uncompress the file with 
```
tar zxvf rascas-master.tar.gz
```
And rename the directory (to be consistent with the git method)
```
mv rascas-master rascas
```

## Building the code

Now, you can go to `rascas/f90/` and compile the code with
```
cd rascas/f90
make all F90=mpif90
```

If you don't have MPI (or if you have a single core processor), use the following instruction to compile single-processor versions of the codes
```
make all F90=gfortran MPI=0
```

## Testing the installation (optional)

To test your installation, you can run the first tutorial notebook
```
cd rascas/tutorials/IdealisedModels/tutorial_1/
make all F90=mpif90
jupyter-notebook Tutorial-1.ipynb
```

If you don’t have MPI, you should set ``useMPI=False`` in the first cell.

You should be able to run succesfully the whole notebook!


***

**Contact for questions or comments:** rascas@univ-lyon1.fr

***

