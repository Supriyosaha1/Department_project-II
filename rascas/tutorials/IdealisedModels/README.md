This directory contains notebooks to illustrate how to use idealised models to describe the ISM/CGM medium within **RASCAS**. 
Each tutorial implements a specific model (analytic distribution and kinematics of gas) described in the file `module_idealised_model.f90`. These modules can be used as templates to develop your own model. Since each module is different, **RASCAS** has to be recompiled in each directory.

- **tutorial_1:** continuum through a static and homogeneous sphere

- **tutorial_2:** emission line through a static and homogeneous sphere

- **tutorial_3:** model of Burchett (thick shell with density $\propto 1/r^2$  and velocity going linearly with radius).

- **tutorial_4:** model of Burchett with cold streams ... (uses mock cubes)

These notebooks served for the hands-on session during the [Saas Fee school 2023](https://www.astro.unige.ch/saasfee2023/).  
