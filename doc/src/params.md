
# Runtime parameters for rascas

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




#### `[LyaPhotonsFromGas]`
The section `[LyaPhotonsFromGas]` is used to configure the generation of Lyman-alpha photon packets from gas emission processes.

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


#### `[mesh]`
The section `[mesh]` is used to define mesh refinement parameters and settings.
| Parameter              | Default Value | Fortran type | Description        |
|:-----------------------|:--------------|:-------------|:-------------------|
| `verbose`              | `T`           | `logical`    | Set verbosity to True or False for mesh operations |
| `refine_lmax`          | `8`           | `integer`    | Maximum refinement level for adaptive mesh refinement |
| `refine_err_grad_d`    | `-1.0d0`      | `real`       | Error threshold for density gradient-based refinement |
| `refine_err_grad_v`    | `-1.d0`       | `real`       | Error threshold for velocity gradient-based refinement |
| `refine_dv_over_vth`   | `F`           | `logical`    | Enable/disable velocity difference threshold relative to thermal velocity for refinement |


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
