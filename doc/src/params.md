
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



