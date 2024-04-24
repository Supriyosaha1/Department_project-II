These tutorials show how to use RASCAS to process an output of a RHD RAMSES simulation in order to generate spectra, images, and data-cubes around the H$\alpha$ line. They model both the line emission and the stellar continuum, and accounts for light propagation through the dusty ISM. 

In practise, they show step-by-step how to prepare the parameter files, run the RASCAS suite, and visualize the results. 

The RASCAS suite works in three steps: 
 - Generate initial conditions for the sources of radiation, i.e. generate photon packets which sample each type of source (star particles or gas).  
 - Extract from RAMSES outputs the gas (and dust) distribution in the computational volume. 
 - Actually compute the radiation transfer.
