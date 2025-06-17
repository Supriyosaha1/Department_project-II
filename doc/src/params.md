
# Runtime parameters for rascas

Configuration files are organised in sections that each contains parameters related to a given code module. At runtime, codes from the RASCAS suite will search for the sections they need in the configuration file and parse parameter names. A unique configuration file can thus be used for several codes, and sections which are not needed will be ignored. Sections do not need to be in any particular order.  

Note on paths: we recommend defining absolute paths in the parameters. 

#### `[RASCAS]`
The section `[RASCAS]` is used to define runtime parameters of the rascas code.


| Parameter      | Default Value       | Fortran type | Description        |
|:---------------|:--------------------|:-------------|:-------------------|
| `DomDumpDir`   | `test/ `            | `character`  | Path to the directory where the domain and mesh files are |
| `PhotonICfile` | `Photon_IC_file.dat`| `character`  | Path to the photon packet initial conditions file|
| `fileout`      | `photons_done.dat`  | `character`  | Path to the standard output file|
| `nbundle`      | `10`                | `integer`    | Number of photon packets sent by the master to each worker at each message (use a large value (e.g. 10,000) when there are few scatterings per photon packet, and a low value (e.g. 10) for Lyman-alpha) | 
| `verbose`      | `T`                 | `logical`    | set verbosity to True or False |



#### `[gas_composition]`
The section `[gas_composition]` is used to define the collection of absorbers. It is read in `module_gas_composition.f90`.


| Parameter | Default Value | Fortran type | Description     |
|:-------------------|:----------------------|:------------|:------------------- |
| `nscatterer`       | `1`                   | `integer`   | Number of scatterers in the run    |
| `scatterer_names`  | `HI-1216`             | `character` | List of names of scatterers (e.g., lines)    |
| `atomic_data_dir`  | `../ions_parameters/` | `character` | Directory where the atomic data files are located     |
| `krome_data_dir`   | `./'`                 | `character` | Directory where Krome metallic ion density files are located    |
| `f_ion`            | `0.01`                | `real`      | Ionization factor used in dust density computation    |
| `Zref`             | `0.005`               | `real`      | Reference metallicity (e.g., ~0.005 for SMC, ~0.01 for LMC)|
| `vturb_kms`        | `20.0`                | `real`      | Constant turbulent velocity across the simulation, in km/s          |
| `ignoreDust`       | `.false.`             | `logical`   | If true, dust is ignored during radiative transfer (but still defined)    |


#### `[IdealisedModel]`
The section `[IdealisedModel]` contains parameters that define the idealized models. This section is read by `module_idealised_model.f90`.

| Parameter                      | Default Value | Fortran type  | Description |
|:-------------------------------|:--------------|:--------------|:------------|
| `ColumnDensity_cgs`            | `1d15`        | `real`       | Column density from the center to the edge of the sphere [cm^-2] |
| `idealised_model_box_size_cm`  | `1d24`        | `real`       | Physical size of the simulation box [cm] |
| `Radius_boxUnits`              | `0.48`        | `real`       | Radius of the sphere in box-size units |
| `Temperature`                  | `1d4`         | `real`       | Gas temperature [K] |
| `TurbulentVelocity_kms`        | `10.0`        | `real`       | Turbulent velocity dispersion [km/s] |



#### `[master]`
The section `[master]` is used to define global runtime parameters and checkpoint/restart options.

| Parameter         | Default Value          | Fortran type      | Description |
|:------------------|:-----------------------|:------------------|:------------|
| `verbose`         | `.false.`              | `logical`         | Set verbosity flag |
| `restart`         | `.false.`              | `logical`         | If True, resume run from backup file specified by `PhotonBakFile` |
| `PhotonBakFile`   | `backup_photons.dat`   | `character`       | Path to the photon backup file |
| `dt_backup`       | `7200.`                | `real`            | Time interval between backups [seconds] |



#### `[mesh]`
The section `[mesh]` is used to define mesh refinement parameters and settings.
| Parameter              | Default Value | Fortran type | Description        |
|:-----------------------|:--------------|:-------------|:-------------------|
| `verbose`              | `T`           | `logical`    | Set verbosity to True or False for mesh operations |
| `refine_lmax`          | `8`           | `integer`    | Maximum refinement level for adaptive mesh refinement |
| `refine_err_grad_d`    | `-1.0d0`      | `real`       | Error threshold for density gradient-based refinement |
| `refine_err_grad_v`    | `-1.d0`       | `real`       | Error threshold for velocity gradient-based refinement |
| `refine_dv_over_vth`   | `F`           | `logical`    | Enable/disable velocity difference threshold relative to thermal velocity for refinement |



#### `[mock]`
The section `[mock]` is used to define parameters for the generation of mock observations.

| Parameter             | Default Value | Fortran type      | Description |
|:----------------------|:--------------|:------------------|:------------|
| `nDirections`         | `0`           | `integer`         | Number of viewing directions for mock observations |
| `mock_parameter_file` | —            | `character`   | Path to the mock parameter configuration file |
| `mock_outputfilename` | —            | `character`   | Prefix for output mock files (including absolute path); suffixes `.flux`, `.image`, `.spectrum`, or `.cube` will be added automatically |

The mock parameter configuration file is a text file that contains one block per direction (or *pointing*). Each block must contain one line per entry as detailed in the table below. Each block of this file is read by the routine `read_a_mock_param_set` from `module_mock.f90`. 

| Line number | description | units | example | comment |
|:------------|:------------|:------|:--------|:--------|
|1 | a normalised vector giving the direction of observation | none |  `0., 0., 1.` | |
|2 | the coordinates of the target | code units |  `0.5, 0.243, 0.12`| | 
|3 | the radius of the circular aperture within which we collect flux | code units | `0.01` | If the value is 0, no flux is computed. A non-zero value will generate a `.flux` output. | 
|4 | number of pixels, aperture radius, min and max wavelengths | no units, code units, Angstrom, Angstrom | `200, 0.01, 1300., 1600.` | If the number of pixels is set to zero, no spectrum is computed. Otherwise, a spectrum is computed and written in the `.spectrum` file.|   
|5 | number of pixels on each side of the image, size of the image   | no units, code units | `200, 0.01`| If the number of pixels is zero, no image is computed. Otherwise, an image is computed and written to the `.image` file. | 
|6 | number of pixels in the spectral dimension, number of pixels in the spatial dimensions, min and max wavelengths, side | none, none, Angstrom, Angstrom, code units | `50, 100, 1200, 1230, 0.01`| If any number of pixels is zero, no cube is computed. Otherwise, a cube is computed and written to the `.cube` file. |  



#### `[ramses]`
The section `[ramses]` is used to define RAMSES simulation parameters and RT variable settings.
| Parameter                           | Default Value | Fortran type | Description        |
|:------------------------------------|:--------------|:-------------|:-------------------|
| `self_shielding`                    | `T`           | `logical`    | If true, reproduce self-shielding approximation made in RAMSES to compute nHI |
| `ramses_rt`                         | `F`           | `logical`    | If true, read RAMSES-RT output and compute nHI and T accordingly |
| `read_rt_variables`                 | `F`           | `logical`    | If true, read RT variables (e.g. to compute heating terms) |
| `cosmo`                             | `T`           | `logical`    | If false, assume idealised simulation |
| `use_proper_time`                   | `F`           | `logical`    | If true, use proper time instead of conformal time for cosmo runs |
| `particle_families`                 | `F`           | `logical`    | If true, all particles have an extra family field, and the header file is different |
| `verbose`                           | `F`           | `logical`    | Display some run-time info on this module |
| `itemp`                             | `5`           | `integer`    | Index of thermal pressure |
| `imetal`                            | `6`           | `integer`    | Index of metallicity |
| `ihii`                              | `7`           | `integer`    | Index of HII fraction |
| `iheii`                             | `8`           | `integer`    | Index of HeII fraction |
| `iheiii`                            | `9`           | `integer`    | Index of HeIII fraction |
| `deut2H_nb_ratio`                   | `3.0d-5`      | `real`       | Ratio between deuterium and hydrogen |
| `recompute_particle_initial_mass`   | `F`           | `logical`    | If true, recompute particle initial mass |
| `tdelay_SN`                         | `10.`         | `real`       | Time delay for supernovae in Myr |
| `recyc_frac`                        | `0.8`         | `real`       | Recycling fraction to correct for mass of stars formed |



#### `[uparallel]`
The section `[uparallel]` is used to define parameters for parallel velocity computation methods.
| Parameter      | Default Value | Fortran type | Description        |
|:---------------|:--------------|:-------------|:-------------------|
| `method`       | `RASCAS`      | `character`  | Computation method, may be 'Smith', 'Semelin', or 'RASCAS' |
| `xForGaussian` | `8.0`         | `real`       | Above this value, use a Gaussian to draw u_parallel |



#### `[RASCAS-serial]`
The section `[RASCAS-serial]` is used to define input/output files and general settings for a serial RASCAS run.

| Parameter        | Default Value           | Fortran type      | Description |
|:------------------|:------------------------|:-----------------|:------------|
| `DomDumpDir`      | `test/`                 | `character`      | Directory where the outputs of `CreateDomDump` are located |
| `PhotonICFile`    | `Photon_IC_file.dat`    | `character`      | Path to the file containing the initial photon packets |
| `fileout`         | `photons_done.dat`      | `character`      | Path to the output file containing processed photons |
| `verbose`         | `.true.`                | `logical`        | Set verbosity flag |



#### `[CreateDomDump]`
The section `[CreateDomDump]` is used to configure the generation of the meshes (or domain dumps) from a simulation output. It is read by the code `CreateDomDump.f90`.

| Parameter        | Default Value           | Fortran type      | Description |
|:------------------|:------------------------|:-----------------|:------------|
| `DomDumpDir`         | `test/`                 | `character`      | Directory where the outputs of `CreateDomDump` will be written |
| `repository`         | `./`                       | `character`      | Ramses run directory (where all output_xxxxx dirs are) |
| `snapnum`               | `1`                         | `integer`          | Ramses output number to use |
| `reading_method` | `fullbox`             | `character`      | strategy to read ramses data, could be `fullbox`, `hilbert`, `select\_onthefly`, or `select\_onthefly\_h` |
| `comput_dom_type`       | `sphere`                   | `character`     | Type of the computational domain (e.g. cube, sphere, shell, slab)|
| `comput_dom_pos`         | `(0.5, 0.5, 0.5)` | `real`          | Center position of the computational domain [code units] |
| `comput_dom_rsp`         | `0.3`                         | `real`          | Radius of the spherical computational domain [code units] |
| `comput_dom_size`       | `0.3`                         | `real`          | Size of the cubic computational domain [code units] |
| `comput_dom_rin`         | `0.0`                         | `real`          | Inner radius of a shell computational domain [code units] |
| `comput_dom_rout`       | `0.3`                         | `real`          | Outer radius of a shell computational domain [code units] |
| `comput_dom_thickness`  | `0.1`                    | `real`          | Thickness of a slab computational domain [code units] |
| `decomp_dom_type`       | `sphere`              | `character` | Type of the decomposition domain(s) (e.g. cube, sphere, shell, slab). Only works with a number of domains, sharing the same type. |
| `decomp_dom_ndomain` | `1`                        | `integer`    | Number of domains in the decomposition| 
| `decomp_dom_xc`           | `0.5`                    | `real`          | x center of domain(s) [code units]. A list of `decomp_dom_ndomain` positions should be provided, e.g. `0.25, 0.5, 0.75`| 
| `decomp_dom_yc`           | `0.5`                    | `real`          | y center of domain(s) [code units] |
| `decomp_dom_zc`           | `0.5`                    | `real`          | z center of domain(s) [code units] |
| `decomp_dom_rsp`         | `0.3`                    | `real`          | Radius of the spherical domain(s) [code units] |
| `decomp_dom_size`       | `0.3`                    | `real`          | Size of the cubic domain(s) [code units] |
| `decomp_dom_rin`         | `0.0`                    | `real`          | Inner radius of the shell domain(s) [code units] |
| `decomp_dom_rout`       | `0.3`                    | `real`          | Outer radius of the shell domain(s) [code units] |
| `decomp_dom_thickness`  | `0.1`               | `real`          | Thickness of the slab domain(s) [code units] |
| `verbose`         | `.true.`                | `logical`        | Set verbosity flag |




#### `[LyaPhotonsFromGas]`
The section `[LyaPhotonsFromGas]` is used to configure the generation of Lyman-alpha photon packets from gas emission processes. It is read by the code `LyaPhotonsFromGas.f90`.

| Parameter | Default Value | Fortran type | Description |
|:----------|:--------------|:-------------|:------------|
| `outputfileRec` | `OutputPath/IC_recomb.dat`                 | `character` | Path to the output file for recombination photon packets |
| `outputfileCol` | `OutputPath/IC_coll.dat`                 | `character` | Path to the output file for collisional photon packets |
| `repository` |`simulation_path/`                 | `character` | Path to the base repository directory containing simulation data |
| `snapnum` | `099`                 | `integer` | Snapshot number to process from the simulation |
| `emission_dom_type` |`cube`                 | `character` | Type of emission domain geometry (e.g., cube, sphere) |
| `emission_dom_pos` |`5e-01 5e-01 5e-01`                 | `real(3)` | Position coordinates (x, y, z) of the emission domain center |
| `emission_dom_rsp` |`1e-02`                 | `real` | Spherical domain radius parameter |
| `emission_dom_rin` |`1e-03`                 | `real` | Inner radius for shell-type emission domains |
| `emission_dom_rout` |`1e-02`                 | `real` | Outer radius for shell-type emission domains |
| `emission_dom_size` |`1e-02`                 | `real` | Size parameter of the emission domain |
| `emission_dom_thickness` |`1e-02`                 | `real` | Thickness parameter for shell-type emission domains |
| `nPhotonPackets` | `1000000`                 | `integer` | Total number of photon packets to generate |
| `ranseed` |`-100`                 | `integer` | Random seed for photon packet generation |
| `doRecombs` | `T`                 | `logical` | Enable processing of recombination photons |
| `doColls` | `F`                | `logical` | Enable processing of collisional photons |
| `tcool_resolution` | `5.0`                 | `real` | Temperature resolution factor for cooling calculations |
| `verbose` | `T` | `logical`                 | Set verbosity to True or False |



#### `[HaPhotonsFromGas]`
The section `[HaPhotonsFromGas]` is used to define runtime parameters related to the generation of H-alpha photons from gas emission processes. It is read by the code `HaPhotonsFromGas.f90`. 

| Parameter                 | Default Value         | Fortran type    | Description |
|:--------------------------|:----------------------|:----------------|:------------|
| `outputfileRec`           | `HaPhotIC.rec`        | `character`     | Path to file to which recombination photons will be written |
| `outputfileCol`           | `HaPhotIC.col`        | `character`     | Path to file to which collisional photons will be written |
| `repository`              | `./`                  | `character`     | RAMSES run directory (where all `output_xxxxx` folders are) |
| `snapnum`                 | `1`                   | `integer`       | RAMSES snapshot number to use |
| `emission_dom_type`       | `sphere`              | `character`     | Type of the emission domain (e.g. cube, sphere, shell, slab)|
| `emission_dom_pos`        | `(0.5, 0.5, 0.5)`     | `real`          | Center position of the emission domain [code units] |
| `emission_dom_rsp`        | `0.3`                 | `real`          | Radius of the spherical emission domain [code units] |
| `emission_dom_size`       | `0.3`                 | `real`          | Size of the cubic emission domain [code units] |
| `emission_dom_rin`        | `0.0`                 | `real`          | Inner radius of a shell emission domain [code units] |
| `emission_dom_rout`       | `0.3`                 | `real`          | Outer radius of a shell emission domain [code units] |
| `emission_dom_thickness`  | `0.1`                 | `real`          | Thickness of a slab emission domain [code units] |
| `nphotons`                | `10000`               | `integer`       | Number of photon packets to generate |
| `ranseed`                 | `-100`                | `integer`       | Random number generator seed for photon sampling |
| `doRecombs`               | `.false.`             | `logical`       | Enable sampling of emission from recombinations |
| `doColls`                 | `.false.`             | `logical`       | Enable sampling of emission from collisional excitations |
| `tcool_resolution`        | `0.0`                 | `real`          | Cells with `dt * tcool_resolution > tcool` are excluded from collisional emission |
| `verbose`                 | `.true.`              | `logical`       | Set verbosity flag |



#### `[GenerateAMRmodel]`
The section `[GenerateAMRmodel]` is used to configure the generation of the meshes (or domain dumps) for idealised (analytical) models. It is read by the code `GenerateAMRmodel.f90`.

| Parameter        | Default Value           | Fortran type      | Description |
|:------------------|:------------------------|:-----------------|:------------|
| `DomDumpDir`         | `test/`                 | `character`      | Directory where the outputs of `GenerateAMRmodel` will be written |
| `comput_dom_type`       | `sphere`                   | `character`     | Type of the computational domain (e.g. cube, sphere, shell, slab)|
| `comput_dom_pos`         | `(0.5, 0.5, 0.5)` | `real`          | Center position of the computational domain [code units] |
| `comput_dom_rsp`         | `0.3`                         | `real`          | Radius of the spherical computational domain [code units] |
| `comput_dom_size`       | `0.3`                         | `real`          | Size of the cubic computational domain [code units] |
| `comput_dom_rin`         | `0.0`                         | `real`          | Inner radius of a shell computational domain [code units] |
| `comput_dom_rout`       | `0.3`                         | `real`          | Outer radius of a shell computational domain [code units] |
| `comput_dom_thickness`  | `0.1`                    | `real`          | Thickness of a slab computational domain [code units] |
| `decomp_dom_type`       | `sphere`              | `character` | Type of the decomposition domain(s) (e.g. cube, sphere, shell, slab). Only works with a number of domains, sharing the same type. |
| `decomp_dom_ndomain` | `1`                        | `integer`    | Number of domains in the decomposition| 
| `decomp_dom_xc`           | `0.5`                    | `real`          | x center of domain(s) [code units]. A list of `decomp_dom_ndomain` positions should be provided, e.g. `0.25, 0.5, 0.75`| 
| `decomp_dom_yc`           | `0.5`                    | `real`          | y center of domain(s) [code units] |
| `decomp_dom_zc`           | `0.5`                    | `real`          | z center of domain(s) [code units] |
| `decomp_dom_rsp`         | `0.3`                    | `real`          | Radius of the spherical domain(s) [code units] |
| `decomp_dom_size`       | `0.3`                    | `real`          | Size of the cubic domain(s) [code units] |
| `decomp_dom_rin`         | `0.0`                    | `real`          | Inner radius of the shell domain(s) [code units] |
| `decomp_dom_rout`       | `0.3`                    | `real`          | Outer radius of the shell domain(s) [code units] |
| `decomp_dom_thickness`  | `0.1`               | `real`          | Thickness of the slab domain(s) [code units] |
| `verbose`         | `.true.`                | `logical`        | Set verbosity flag |

