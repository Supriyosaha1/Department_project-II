
# Runtime parameters for rascas

The section `[RASCAS]` is used to define runtime parameters of the main rascas code.


| Variable name, syntax, default value | Fortran type  | Description               |
|:---------------------------- |:------------- |:------------------------- |
| `DomDumpDir = test/ `        |  `character`  | Directory where the domain and mesh files are |
| `PhotonICfile = Photon_IC_file.dat` | `character` | name of the Photon Packet Initial Conditions file (with path) |
| `fileout = photons_done.dat`  | `character` | name of the output file |
| `nbundle = 10`               | `Integer`   | number of PP sent by the master to each worker at each message | 
| `verbose = T`                | `Logical`   | set verbosity to True of False|


Alternatively, here is the ``[RASCAS]`` section to be included in the parameter file.


```
[RASCAS]
  DomDumpDir   = ./test/               # directory where the domain and mesh files are
  PhotonICFile = ppic.dat              # name of the Photon Packet Initial Conditions file (with path)
  fileout      = photons_done.dat      # name of the output file
  nbundle      = 10                    # 
  verbose      = T                     # 
```




