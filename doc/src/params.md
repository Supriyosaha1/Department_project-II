
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
| `emission_dom_pos` |`1e-01 1e-01 1e-01`                 | `real(3)` | Position coordinates (x, y, z) of the emission domain center |
| `emission_dom_rsp` |`5e-03`                 | `real` | Spherical domain radius parameter |
| `emission_dom_rin` |`2e-03`                 | `real` | Inner radius for shell-type emission domains |
| `emission_dom_rout` |`4e-03`                 | `real` | Outer radius for shell-type emission domains |
| `emission_dom_size` |`1e-02`                 | `real` | Size parameter of the emission domain |
| `emission_dom_thickness` |`1e-02`                 | `real` | Thickness parameter for shell-type emission domains |
| `nPhotonPackets` | `1000000`                 | `integer` | Total number of photon packets to generate |
| `ranseed` |`-100`                 | `integer` | Random seed for photon packet generation |
| `doRecombs` | `T`                 | `logical` | Enable processing of recombination photons |
| `doColls` | `F`                | `logical` | Enable processing of collisional photons |
| `tcool_resolution` | `5.0`                 | `real` | Temperature resolution factor for cooling calculations |
| `verbose` | `T` | `logical`                 | Set verbosity to True or False |
