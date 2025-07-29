
# Runtime parameters for rascas

Configuration files are organised in sections that each contains parameters related to a given code module. At runtime, codes from the RASCAS suite will search for the sections they need in the configuration file and parse parameter names. A unique configuration file can thus be used for several codes, and sections which are not needed will be ignored. Sections do not need to be in any particular order. Examples of parameter files can be found in the tests directory. 

Note on paths: we recommend defining absolute paths in the parameters. 

#### `[RASCAS]`
The section `[RASCAS]` is used to define runtime parameters of the rascas code.


| Parameter      | Default Value       | Fortran type | Description        |
|:---------------|:--------------------|:-------------|:-------------------|
| `DomDumpDir`   | `test/ `            | `character`  | Path to the directory where the domain and mesh files are |
| `PhotonICfile` | `Photon_IC_file.dat`| `character`  | Path to the photon packet initial conditions file|
| `fileout`      | `photons_done.dat`  | `character`  | Path to the standard output file|
| `nbundle`      | `10`                | `integer`    | Number of photon packets sent by the master to each worker at each message (use a large value (e.g. 10,000) when there are few scatterings per photon packet, and a low value (e.g. 10) for Lyman-alpha) | 
| `verbose`      | `.true.`                 | `logical`    | set verbosity to True or False |



#### `[dust]`
The section `[dust]` contains parameters that define the dust modelling. This section is read by `module_dust_model.f90`. 

| Parameter             | Default Value | Fortran type   | Description |
|:--------------|:-----------|:------------|:------------|
| `albedo`            | `0.32`          | `real`           | dust albedo, the value is dust and wavelength dependant. Default 0.32 for Lya, from Li & Draine 2001 | 
| `g_dust`            | `0.73`          | `real`           | g parameter of the Henyey-Greenstein phase function for dust scattering. Default 0.73 from Li & Draine 2001 |
| `dust_model`    | `SMC`            | `character` | Dust extinction law. Could be SMC or LMC |



#### `[gas_composition]`
The section `[gas_composition]` is used to define the collection of absorbers that will interact with photons. It is read in `module_gas_composition.f90` and uses atomic data from files in the directory `atomic_data_dir` which are described below. 

| Parameter | Default Value | Fortran type | Description     |
|:-------------------|:----------------------|:------------|:------------------- |
| `nscatterer`       | --                   | `integer`   | Number of scatterers in the run    |
| `scatterer_names`  | --           | `character` | List of names of scatterers (e.g., lines)    |
| `atomic_data_dir`  | `../ions_parameters/` | `character` | Directory where the atomic data files are located     |
| `krome_data_dir`   | `./'`                 | `character` | Directory where Krome metallic ion density files are located    |
| `f_ion`            | `0.01`                | `real`      | Ionization factor used in dust density computation    |
| `Zref`             | `0.005`               | `real`      | Reference metallicity (e.g., ~0.005 for SMC, ~0.01 for LMC)|
| `vturb_kms`        | `20.0`                | `real`      | Constant turbulent velocity across the simulation, in km/s          |
| `ignoreDust`       | `.false.`             | `logical`   | If true, dust is ignored during radiative transfer (but still defined)    |

The photons may interact with dust and scatterers. Scatterers are defined by their absorption chanel (e.g. `HI-1216` refers to the absorption of radiation at H Lyman-alpha). Some absorption channels may lead to multiple de-excitation routes (e.g. fluorescent emission from SiII). This is allowed and described in the atomic data files. There are example files in the `ions_parameters` directory shipped with the code. The generic format is as follows: 
```
! SiII - 1260.42 line
m_ion          = 28.085d0.      ! mass in atomic mass unit
name_ion       = SiII           ! name of the ion (prefix of line name)
lambda_cm      = 1260.42d-8     ! wavelength in cm
A              = 2.57d9         ! Einstein coefficient A21
f              = 1.22d0.        ! Oscillator strength
n_fluo         = 1              ! number of fluorescent channels (i.e. secondary de-escitation channels)
lambda_fluo_cm = 1265.02d-8     ! wavelengths of these channels 
A_fluo         = 4.73d8         ! Einseint coeff of these channels
```

HI-Lyman-alpha is a specific case which has the extra following parameters: 
```
isotropic = F        ! controls whether scattering is isotropic or not. 
recoil    = T        ! controls whether the recoil effect is included or not
core_skip = F        ! allows the core-skipping approximation
xcritmax  = 100      ! parameter to moderate the core-skipping approximation
lee_correction = F   ! apply the Lee correction. 
```



#### `[IdealisedModel]`
The section `[IdealisedModel]` contains parameters that define the idealized models. This section is read by `module_idealised_model.f90`.

| Parameter                      | Default Value | Fortran type  | Description |
|:-------------------------------|:--------------|:--------------|:------------|
| `ColumnDensity_cgs`            | `1d15`        | `real`       | Column density from the center to the edge of the sphere [cm^-2] |
| `box_size_cm`  | `1d24`        | `real`       | Physical size of the simulation box [cm] |
| `Radius_boxUnits`              | `0.48`        | `real`       | Radius of the sphere in box-size units |
| `Temperature`                  | `1d4`         | `real`       | Gas temperature [K] |
| `TurbulentVelocity_kms`        | `10.0`        | `real`       | Turbulent velocity dispersion [km/s] |



#### `[master]`
The section `[master]` is used to define global runtime parameters and checkpoint/restart options. This section is read by `module_master.f90`.

| Parameter         | Default Value          | Fortran type      | Description |
|:------------------|:-----------------------|:------------------|:------------|
| `verbose`         | `.false.`              | `logical`         | Set verbosity flag |
| `restart`         | `.false.`              | `logical`         | If True, resume run from backup file specified by `PhotonBakFile` |
| `PhotonBakFile`   | `backup_photons.dat`   | `character`       | Path to the photon backup file |
| `dt_backup`       | `7200.`                | `real`            | Time interval between backups [seconds] |



#### `[mesh]`
The section `[mesh]` is used to define mesh refinement parameters and settings.  It is read in `module_mesh.f90`. The parameters `refine_lmax`, `refine_err_grad_d`, `refine_err_grad_v`, and `refine_dv_over_vth` are only used to generate AMR grids for idealised models (i.e. with `GenerateAMRmodel`). These four parameters are ignored when reading simulation outputs. 
| Parameter                              | Default Value       | Fortran type       | Description        |
|:-----------------------|:--------------|:-------------|:-------------------|
| `verbose`                            | `.true.`            | `logical`        | Set verbosity to True or False for mesh operations |
| `refine_lmax`                    | `8`                      | `integer`        | Maximum level of refinement allowed in the mesh.|
| `refine_err_grad_d`        | `-1.0d0`            | `real`              | Parameter to control the refinement on density gradients. A cell will be refined if {math}`(\rho_{max}-\rho_{min})/(\rho_{max}+\rho_{min}) >`  `refine_err_grad_d`, where {math}`\rho_{max}` and {math}`\rho_{min}` are the max and min values of densities found accross the cell. Setting `refine_err_grad_d=0.2` triggers refinement when the density varies by more than 10% accross a cell. If a negative value is passed, this criterion is ignored. |
| `refine_err_grad_v`        | `-1.d0`              | `real`              | Parameter to control the refinement on velocity gradients. Same as for density above, but now computed on the norm of the velocity field. If a negative value is passed, this criterion is ignored.  |
| `refine_dv_over_vth`      | `.false.`          | `logical`        | Parameter to control the refinement on the velocity field. If set to True, a cell is refined when variations of the norm of the velocity accross the cell are larger than the local thermal velocity. |



#### `[mock]`
The section `[mock]` is used to define parameters for the generation of mock observations. This section is read by `module_mock.f90`.

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
The section `[ramses]` is used to define RAMSES simulation parameters and RT variable settings. This section is read by `module_ramses.f90`.
| Parameter                           | Default Value      | Fortran type   | Description        |
|:------------------------------------|:-----------|:-------------|:-------------------|
| `self_shielding`           | `.true.`           | `logical`    | If true, reproduce self-shielding approximation made in RAMSES to compute nHI |
| `ramses_rt`                     | `.false.`         | `logical`    | If true, read RAMSES-RT output and compute nHI and T accordingly |
| `read_rt_variables`     | `.false.`         | `logical`    | If true, read RT variables (e.g. to compute heating terms) |
| `cosmo`                             | `.true.`           | `logical`    | If false, assume idealised simulation |
| `use_proper_time`         | `.false.`         | `logical`    | If true, use proper time instead of conformal time for cosmo runs |
| `particle_families`     | `.false.`         | `logical`    | If true, all particles have an extra family field, and the header file is different |
| `verbose`                         | `.false.`         | `logical`    | Display some run-time info on this module |
| `itemp`                             | `5`           | `integer`    | Index of thermal pressure |
| `imetal`                           | `6`           | `integer`    | Index of metallicity |
| `ihii`                               | `7`           | `integer`    | Index of HII fraction |
| `iheii`                             | `8`           | `integer`    | Index of HeII fraction |
| `iheiii`                           | `9`           | `integer`    | Index of HeIII fraction |
| `deut2H_nb_ratio`         | `3.0d-5`      | `real`       | Ratio between deuterium and hydrogen |
| `recompute_particle_initial_mass`   | `.false.` | `logical`    | If true, recompute particle initial mass |
| `tdelay_SN`                     | `10.`         | `real`       | Time delay for supernovae in Myr |
| `recyc_frac`                   | `0.8`         | `real`       | Recycling fraction to correct for mass of stars formed |



#### `[uparallel]`
The section `[uparallel]` is used to define parameters for parallel velocity computation methods. This section is read by `module_uparallel.f90`.
| Parameter             | Default Value        | Fortran type     | Description        |
|:---------------|:--------------|:------------|:-------------------|
| `method`             | `RASCAS`            | `character`  | Computation method, may be 'Smith', 'Semelin', or 'RASCAS' |
| `xForGaussian` | `8.0`                  | `real`            | Above this value, use a Gaussian to draw u_parallel |



#### `[voigt]`
The section `[voigt]` is used to define the numerical approximation to the Voigt function, used at each interaction. This section is read by `module_voigt.f90`.

| Parameter        | Default Value           | Fortran type      | Description |
|:------------------|:------------------------|:-----------------|:------------|
| `approximation`      | `COLT`                 | `character`      | could be 'COLT', 'Tasitsiomi' or 'Humlicek_w4' |



#### `[worker]`
The section `[worker]` is used to define runtime parameters of the workers. This section is read by `module_worker.f90`. 

| Parameter        | Default Value           | Fortran type       | Description |
|:------------|:----------------|:-------------|:------------|
|`verbose`       | `.false.`                | `logical`        | Set verbosity flag |



#### `[RASCAS-serial]`
The section `[RASCAS-serial]` is used to define input/output files and general settings for a serial RASCAS run.

| Parameter              | Default Value                       | Fortran type         | Description |
|:---------------|:----------------------|:--------------|:------------|
| `DomDumpDir`      | `test/`                              | `character`      | Directory where the outputs of `CreateDomDump` are located |
| `PhotonICFile`  | `Photon_IC_file.dat`    | `character`      | Path to the file containing the initial photon packets |
| `fileout`            | `photons_done.dat`        | `character`      | Path to the output file containing processed photons |
| `verbose`            | `.true.`                            | `logical`          | Set verbosity flag |



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
| `decomp_dom_rsp`         | `0.35`                  | `real`          | Radius of the spherical domain(s) [code units] |
| `decomp_dom_size`       | --                          | `real`          | Size of the cubic domain(s) [code units] |
| `decomp_dom_rin`         | --                          | `real`          | Inner radius of the shell domain(s) [code units] |
| `decomp_dom_rout`       | --                          | `real`          | Outer radius of the shell domain(s) [code units] |
| `decomp_dom_thickness`  | --                     | `real`          | Thickness of the slab domain(s) [code units] |
| `verbose`         | `.true.`                | `logical`        | Set verbosity flag |




#### `[LyaPhotonsFromGas]`
The section `[LyaPhotonsFromGas]` is used to configure the generation of Lyman-alpha photon packets from gas emission processes. It is read by the code `LyaPhotonsFromGas.f90`.

| Parameter | Default Value | Fortran type | Description |
|:----------|:--------------|:-------------|:------------|
| `outputfileRec` | `LyaPhotIC.recLya`                 | `character` | Path to the output file for recombination photon packets |
| `outputfileCol` | `LyaPhotIC.colLya`                 | `character` | Path to the output file for collisional photon packets |
| `repository` |`./`                 | `character` | Path to the base repository directory containing simulation data |
| `snapnum` | `1`                 | `integer` | Snapshot number to process from the simulation |
| `emission_dom_type` |`sphere`                 | `character` | Type of emission domain geometry (e.g., cube, sphere) |
| `emission_dom_pos` |`5e-01 5e-01 5e-01`                 | `real(3)` | Position coordinates (x, y, z) of the emission domain center |
| `emission_dom_rsp` |`0.3`                 | `real` | Spherical domain radius parameter |
| `emission_dom_rin` |`0.`                 | `real` | Inner radius for shell-type emission domains |
| `emission_dom_rout` |`0.3`                 | `real` | Outer radius for shell-type emission domains |
| `emission_dom_size` |`0.3`                 | `real` | Size parameter of the emission domain |
| `emission_dom_thickness` |`0.1`                 | `real` | Thickness parameter for shell-type emission domains |
| `nPhotonPackets` | `10000`                 | `integer` | Total number of photon packets to generate |
| `ranseed` |`-100`                 | `integer` | Random seed for photon packet generation |
| `doRecombs` | `.false.`                 | `logical` | Enable processing of recombination photons |
| `doColls` | `.true.`                | `logical` | Enable processing of collisional photons |
| `tcool_resolution` | `3.0`                 | `real` | Temperature resolution factor for cooling calculations |
| `verbose` | `.true.` | `logical`                 | Set verbosity to True or False |



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
| `decomp_dom_type`       | `cube`              | `character` | Type of the decomposition domain(s) (e.g. cube, sphere, shell, slab). Only works with a number of domains, sharing the same type. |
| `decomp_dom_ndomain` | `1`                        | `integer`    | Number of domains in the decomposition| 
| `decomp_dom_xc`           | `0.5`                    | `real`          | x center of domain(s) [code units]. A list of `decomp_dom_ndomain` positions should be provided, e.g. `0.25, 0.5, 0.75`| 
| `decomp_dom_yc`           | `0.5`                    | `real`          | y center of domain(s) [code units] |
| `decomp_dom_zc`           | `0.5`                    | `real`          | z center of domain(s) [code units] |
| `decomp_dom_rsp`         | --                          | `real`          | Radius of the spherical domain(s) [code units] |
| `decomp_dom_size`       | `1.`                      | `real`          | Size of the cubic domain(s) [code units] |
| `decomp_dom_rin`         | --                          | `real`          | Inner radius of the shell domain(s) [code units] |
| `decomp_dom_rout`       | --                          | `real`          | Outer radius of the shell domain(s) [code units] |
| `decomp_dom_thickness`  | --                     | `real`          | Thickness of the slab domain(s) [code units] |
| `verbose`         | `.true.`                | `logical`        | Set verbosity flag |



#### `[PhotonsFromStars]`
The section `[PhotonsFromStars]` is used to configure the generation of continuum photon packets from the stellar populations. It is read by the code `PhotonsFormStars.f90`.

| Parameter    | Default Value | Fortran type | Description |
|:----------|:--------------|:-------------|:------------|
| `outputfile`                     | `PhotICs.dat`             | `character` | Path to the output file for the photon packets |
| `repository`                     |`./`                                | `character` | Path to the base repository directory containing simulation data |
| `snapnum`                           | `1`                                 | `integer`     | Snapshot number to process from the simulation |
| `star_dom_type`               |`cube`                            | `character` | Geometry of the domain containing the emitting stars. Could be cube, sphere, shell, or slab |
| `star_dom_pos`                 |`5e-01 5e-01 5e-01`  | `real(3)`     | Position coordinates (x, y, z) of the domain center |
| `star_dom_rsp`                 |`0.3`                              | `real`           | Spherical domain radius parameter |
| `star_dom_rin`                 |`0.`                                | `real`           | Inner radius for shell-type emission domains |
| `star_dom_rout`               |`0.3`                              | `real`           | Outer radius for shell-type emission domains |
| `star_dom_size`               |`0.3`                              | `real`           | Size parameter of the emission domain |
| `star_dom_thickness`     |`0.1`                              | `real`           | Thickness parameter for shell-type emission domains |
| `spec_type`                       |`Mono`                            | `character` | Spectral shape of the emitting stars. Could be Mono, Gauss, PowLaw, or Table | 
| `spec_SSPdir`                   |`../libs/SSPlibs/`    | `character` | The SSP lib directory |
| `spec_mono_lambda0`       |`1216.`                          | `real`           | Monochromatic emission wavelength [A] |
| `spec_gauss_lambda0`     |`1216.`                          | `real`           | Central wavelength [A] of the Gaussian distribution | 
| `spec_gauss_sigma_kms` |`10.0`                            | `real`           | Line width in velocity [km/s] of the Gaussian distribution | 
| `spec_table_lmin_Ang`   |`1120.`                          | `real`           | min wavelength [A] for the sampling of the tabulated spectrum |
| `spec_table_lmax_Ang`   |`1320.`                          | `real`           | max wavelength [A] for the sampling of the tabulated spectrum | 
| `spec_powlaw_lmin_Ang`                               | `1120.`                          | `real`           | min wavelength [A] for the sampling of a power-law continuum |
| `spec_powlaw_lmax_Ang`                               | `1320.`                          | `real`           | max wavelength [A] for the sampling of a power-law continuum | 
| `spec_powlaw_Fitlmin_Ang`                         | `1100.`                          | `real`           | min wavelength [A] for the fit  |
| `spec_powlaw_Fitlmax_Ang`                         | `1400.`                          | `real`           | max wavelength [A] for the fit | 
| `spec_powlaw _AbsorptionLineClipping` | `.true.`                        | `logical`     | remove absorption lines from fit |
| `nPhotonPackets`             | `1000000`          | `integer`       | Total number of photon packets to generate |
| `ranseed`                           |`-100`                 | `integer`       | Random seed for photon packet generation |
| `verbose`                           | `.true.`            | `logical`       | Set verbosity to True or False |



#### `[ExtractSubvol]`
The code `ExtractSubvol.f90` is a post-processing tool for RAMSES simulations. It is inherited from the code `CreateDomDump.f90`. The basic idea is to extract from a RAMSES snapshot a subvolume of the simulation using the domain strategy of RASCAS. This cutout is then easy to manipulate with Python. It also uses the rebuilding of the oct-tree structure for the selected cells, allowing the use of neighbors to compute some properties (e.g. velocity dispersion). 

| Parameter        | Default Value           | Fortran type      | Description |
|:------------------|:------------------------|:-----------------|:------------|
| `DomDumpDir`         | `test/`                 | `character`      | Directory where the outputs of `CreateDomDump` will be written |
| `repository`         | `./`                       | `character`      | Ramses run directory (where all output_xxxxx dirs are) |
| `snapnum`               | `1`                         | `integer`          | Ramses output number to use |
| `reading_method` | `fullbox`             | `character`      | strategy to read ramses data, could be `fullbox`, `hilbert`, `select_onthefly`, or `select_onthefly_h` |
| `decomp_dom_type`       | `sphere`              | `character` | Type of the decomposition domain (e.g. cube, sphere, shell, slab). |
| `decomp_dom_xc`           | `0.5`                    | `real`          | x center of domain [code units] |
| `decomp_dom_yc`           | `0.5`                    | `real`          | y center of domain [code units] |
| `decomp_dom_zc`           | `0.5`                    | `real`          | z center of domain [code units] |
| `decomp_dom_rsp`         | `0.5`                    | `real`          | Radius of the spherical domain [code units] |
| `decomp_dom_size`       | `0.3`                    | `real`          | Size of the cubic domain [code units] |
| `decomp_dom_rin`         | `0.0`                    | `real`          | Inner radius of the shell domain [code units] |
| `decomp_dom_rout`       | `0.3`                    | `real`          | Outer radius of the shell domain [code units] |
| `decomp_dom_thickness`  | `0.1`               | `real`          | Thickness of the slab domain [code units] |
| `verbose`         | `.false.`                | `logical`        | Set verbosity flag |
| `add_stars`     | `.false.`                | `logical`        | If true add the stellar particles contained in the domain to the cutout. |
| `add_dm`           | `.false.`                | `logical`        | Same for DM particles. Not implemented yet. |


#### `[PhotonsFromSourceModel]`
The section `[PhotonsFromSourceModel]` is used to configure the generation of photon packets from idealised source models. It is read by the code `PhotonsFromSourceModel.f90`.

| Parameter                                  | Default Value                | Fortran type      | Description |
|:-------------------------|:------------------|:-------------|:------------|
| `outputfile`                          | `ppic.dat`                 | `character`   | define the path to the output file |
| `source_type`                        | `pointlike`               | `character`   | can only be point like |
| `source_pos`                          | `0.5 0.5 0.5`           | `real`             | coordinates of the source (x,y,z), in box units  |
| `source_vel`                          | `0.0 0.0 0.0`           | `real`             | velocity of the source (vx,vy,vz), in cm/s    |
| `nphotons`                              | `1000`                         | `integer`       | number of photon to launch from the source |
| `spec_type`                            | `Gauss`                       | `character`   | could be 'Mono' or 'Gauss' or 'PowLaw' or 'Table' or 'TablePowerLaw' |
| `spec_mono_l0_Ang`              | `1215.67`                   | `real`             | wavelength of the source in Angstrom, only used if spec_type=Mono   |
| `spec_gauss_l0_Ang`            | `1215.67`                   | `real`             | central wavelength of the gaussian in Angstrom, only used if spec_type=Gauss   |
| `spec_gauss_sigma_kms`      | `10.0`                         | `real`             | width of the gaussian in km/s, only used if spec_type=Gauss   |
| `spec_powlaw_lmin_Ang`      | `1120.0`                     | `real`             | min wavelength to sample in Angstrom, only used if spec_type=PowLaw   |
| `spec_powlaw_lmax_Ang`      | `1320.0`                     | `real`             | max wavelength to sample in Angstrom, only used if spec_type=PowLaw   |
| `spec_powlaw_beta`              | `-2.3`                         | `real`             | power law index, only used if spec_type=PowLaw   |
| `spec_SSPdir`                        | `../libs/SSPlibs/` | `character`   | path to the SSP lib directory, only used if spec_type=Table or spec_type=TablePowLaw  |
| `spec_table_lmin_Ang`        | `1120.0`                     | `real`             | min wavelength to sample in Angstrom, only used if spec_type=Table   |
| `spec_table_lmax_Ang`        | `1320.0`                     | `real`             | max wavelength to sample in Angstrom, only used if spec_type=Table   |
| `spec_table_age`                  | `10.0`                         | `real`             | age of the stellar population to use in Myr, only used if spec_type=Table or spec_type=TablePowLaw   |
| `spec_table_met`                  | `0.02`                         | `real`             | metallicity of the stellar population to use, only used if spec_type=Table  or spec_type=TablePowLaw  |
| `spec_table_mass`                | `1.e6`                         | `real`             | mass of the source in solar mass, only used if spec_type=Table or spec_type=TablePowLaw   |
| `spec_tpl_lmin_Ang`            | `1120.0`                     | `real`             | min wavelength to sample in Angstrom, only used if spec_type=TablePowLaw   |
| `spec_tpl_lmax_Ang`            | `1320.0`                     | `real`             | max wavelength to sample in Angstrom, only used if spec_type=TablePowLaw   |
| `spec_tpl_Fitlmin_Ang`      | `1100.`                       | `real`             | min wavelength to sample for fit, only used if spec_type=TablePowLaw   |
| `spec_tpl_Fitlmax_Ang`      | `1500.`                       | `real`             | max wavelength to sample for fit, only used if spec_type=TablePowLaw   |
| `spec_tpl_AbsorptionLineClipping` | `.true.`    | `logical`       | if true remove absorption lines from fit, only used if spec_type=TablePowLaw |
| `ranseed`         | `1234`                    | `integer`        | seed for the random generator |
| `verbose`         | `.true.`                | `logical`        | Set verbosity flag |
