# iPWR-SMR-Dynamic-Model

Thermodynamically coupled dynamic model of an iPWR-type small modular reactor (SMR) with a three-region moving-boundary steam generator and physics-based Rankine cycle.

This repository contains the MATLAB, Simulink, and Simscape implementation associated with the integrated dynamic SMR model developed for load-following analysis. The model couples an equation-based NuScale-type integral pressurized water reactor (iPWR), including reactor kinetics, primary-loop thermal hydraulics, and a three-region moving-boundary helical-coil once-through steam generator, with a physics-based secondary Rankine cycle.

The purpose of this repository is to support reproducibility, documentation, and future development of the SMR dynamic modeling framework.

---

## Model Overview

The model consists of two main parts:

### 1. Equation-Based SMR and Steam Generator Model

The MATLAB/Simulink portion represents the nuclear island and steam generator dynamics, including:

* point-kinetics reactor model with one effective delayed neutron precursor group,
* fuel and coolant temperature reactivity feedback,
* lumped fuel and primary-coolant thermal-hydraulic states,
* constant-speed pump-driven primary flow,
* hot-leg and cold-leg transport delays,
* three-region moving-boundary steam generator model:

  * subcooled region, `L1`,
  * saturated/two-phase region, `L2`,
  * superheated region, `L3`,
  * with `L1 + L2 + L3 = LT`.

### 2. Physics-Based Secondary Rankine Cycle

The Simscape/Simscape Fluids portion represents the secondary steam cycle, including:

* steam throttle valve,
* high-pressure and low-pressure turbine stages,
* condenser,
* feedwater pump,
* controlled two-phase reservoir interfaces between the equation-based steam generator and the Simscape network.

This coupling allows the model to capture dynamic back-pressure, pressure-flow interactions, variable turbine enthalpy drop, feedwater/steam flow mismatch, and moving-boundary steam-generator behavior during load-following transients.

## Main Files

### Parameter Definition

`buildParams4.m`

Defines the main physical and numerical parameters used by the reactor, primary loop, and steam generator model. This includes neutron-kinetics parameters, reactivity coefficients, rated thermal power, primary flow rate, steam generator geometry, heat-transfer coefficients, primary-system volumes, material properties, and moving-boundary steam generator constants.

### Steam/Water Property Library

`XSteam.m`

Provides water and steam thermodynamic properties based on the IAPWS IF-97 formulation. This file is used throughout the model for evaluating saturation properties, pressure-temperature properties, pressure-enthalpy properties, density, enthalpy, entropy, internal energy, and heat capacity.

### Steady-State Residual Function

`reactorSteadyResidual2.m`

Defines the nonlinear residual equations used to solve the nominal steady-state operating point of the reactor and three-region steam generator model. This function is called by `steadySolve2.m` through MATLAB's `fsolve`.

### Steady-State Solver

`steadySolve2.m`

Solves the steady-state initialization problem and generates the dynamic initial condition vector. The script saves the results in:

```text
steadyStateResults.mat
```

This file must be generated before running the Simulink/Simscape dynamic model.

### Dynamic ODE/DAE Model

`reactor_ode4.m`

Defines the main dynamic model for the reactor, primary loop, and moving-boundary steam generator. It computes the time derivatives of the 16 dynamic states, including neutron flux, precursor concentration, fuel temperature, coolant temperatures, hot-leg/cold-leg temperatures, steam generator moving-boundary lengths, secondary pressure, outlet enthalpy, tube-wall temperatures, and primary-side steam-generator temperatures.

### Simulink Interface

`reactorBlock.m`

Implements a Level-2 MATLAB S-function that wraps the dynamic MATLAB model for use inside Simulink. It loads the steady-state initial condition file and connects the equation-based SMR/steam-generator model to the Simulink/Simscape secondary cycle.

### Full Simulink/Simscape Model

`reactorBlock_Simulink_Simscape_Full_Controller.slx`

Contains the full integrated dynamic model, including the MATLAB S-function reactor block, Simscape two-phase Rankine-cycle components, and the full valve-pump-rod control architecture.

---

## Required Software

The model was developed using MATLAB/Simulink with Simscape components.

Recommended requirements:

* MATLAB R2023a or newer,
* Simulink,
* Simscape,
* Simscape Fluids,
* Optimization Toolbox.

The Optimization Toolbox is required for `fsolve`, which is used in the steady-state initialization step.

---

## How to Run the Model

### Step 1: Clone the Repository

```bash
git clone https://github.com/USERNAME/iPWR-SMR-Dynamic-Model.git
cd iPWR-SMR-Dynamic-Model
```

Replace `USERNAME` with the appropriate GitHub username or organization name.

### Step 2: Open MATLAB in the Repository Folder

Open MATLAB and set the current folder to the root directory of the repository.

Then add all subfolders to the MATLAB path:

```matlab
clear; clc;
addpath(genpath(pwd));
```

### Step 3: Run the Steady-State Solver

Before running the dynamic Simulink model, solve the nominal steady-state operating point:

```matlab
steadySolve2
```

This script generates:

```text
steadyStateResults.mat
```

The file contains the steady-state solution and the dynamic initial condition vector required by the Simulink S-function.

### Step 4: Open the Full Simulink/Simscape Model

```matlab
open_system('reactorBlock_Simulink_Simscape_Full_Controller.slx');
```

### Step 5: Run the Dynamic Simulation

```matlab
sim('reactorBlock_Simulink_Simscape_Full_Controller.slx');
```

The full model simulates the coupled SMR, moving-boundary steam generator, Rankine cycle, and control system.

---

## Logical Run Order

The files should be used in the following order:

```text
1. XSteam.m
   Dependency for water/steam properties.

2. buildParams4.m
   Builds the parameter structure.

3. reactorSteadyResidual2.m
   Defines the steady-state residual equations.

4. steadySolve2.m
   Solves the nominal steady state and saves steadyStateResults.mat.

5. reactor_ode4.m
   Defines the dynamic reactor and steam-generator model.

6. reactorBlock.m
   Wraps reactor_ode4.m as a Simulink S-function.

7. reactorBlock_Simulink_Simscape_Full_Controller.slx
   Runs the full integrated Simulink/Simscape dynamic model.
```

In practice, the user runs only:

```matlab
steadySolve2
```

followed by the Simulink model.

---

## Dynamic States

The main dynamic model uses 16 states:

```text
1.  phi      normalized neutron flux
2.  C        delayed neutron precursor concentration
3.  T_f      fuel temperature
4.  T_c1     core coolant temperature
5.  T_HL     hot-leg temperature
6.  T_CL     cold-leg temperature
7.  L1       subcooled steam-generator region length
8.  L2       saturated/two-phase steam-generator region length
9.  P_s      secondary-side steam-generator pressure
10. h_out    steam-generator outlet enthalpy
11. T_m3     tube-wall temperature, subcooled region
12. T_m2     tube-wall temperature, saturated region
13. T_m1     tube-wall temperature, superheated region
14. T_p5     primary-side steam-generator temperature, subcooled region
15. T_p3     primary-side steam-generator temperature, saturated region
16. T_p1     primary-side steam-generator temperature, superheated region
```

---

## Initial Conditions and Rated Steady-State Initialization

All dynamic simulations are initialized from the nominal rated operating point of the integrated SMR model. Before running the Simulink/Simscape transient simulation, the MATLAB steady-state solver is used to compute a consistent rated steady-state solution for the reactor, primary loop, and three-region moving-boundary steam generator.

The steady-state solution corresponds to rated operating conditions, including:

* rated reactor thermal power,
* rated primary-loop pressure,
* rated secondary/steam-generator pressure,
* rated primary mass flow rate,
* rated steam mass flow rate,
* nominal hot-leg and cold-leg temperatures,
* nominal steam-generator outlet temperature,
* nominal moving-boundary region lengths `L1`, `L2`, and `L3`.

For the nominal case, the rated steady-state conditions correspond to 160 MW of reactor thermal power, 12.76 MPa of primary pressure, 3.448 MPa of secondary steam pressure, 587 kg/s of primary mass flow rate, and 67.07 kg/s of steam mass flow rate.

The script `steadySolve2.m` solves the nonlinear steady-state residual equations defined in `reactorSteadyResidual2.m`. The resulting solution is saved in:

```text
steadyStateResults.mat
```

This file contains the dynamic initial condition vector `X_dyn0`, which is loaded by `reactorBlock.m` when the Simulink/Simscape model starts. Therefore, the transient load-following simulation begins from a physically consistent rated steady state rather than from manually assigned or arbitrary initial values.

If any major model parameter is changed in `buildParams4.m`, the steady-state solver should be run again before executing the Simulink/Simscape model. This ensures that the initial conditions remain consistent with the updated rated operating point.

---

## Load-Following Case

The main transient case is a 5% reduction in turbine mechanical power demand, from 50 MW to 47.5 MW.

The full controller includes:

* valve control for turbine mechanical power tracking,
* feedwater pump control for steam-generator pressure regulation,
* control-rod reactivity control for primary coolant temperature regulation.

This coordinated valve-pump-rod strategy is used to study the interaction among mechanical power tracking, secondary pressure control, primary temperature regulation, and steam-generator moving-boundary margins.

---

## Notes for Users

1. Run `steadySolve2.m` before running the Simulink model.

2. If physical parameters are changed in `buildParams4.m`, rerun `steadySolve2.m` to regenerate a consistent `steadyStateResults.mat` file.

3. Keep `XSteam.m` on the MATLAB path. The model depends on it for water/steam thermodynamic properties.

4. Be careful with units. The main model uses SI units, while some `XSteam.m` calls use pressure in bar and enthalpy in kJ/kg.

5. The `.slx` file requires Simulink, Simscape, and Simscape Fluids.

---

## Citation

If you use this repository, please cite the associated paper:

```text
Mahboub Rad, A., Jacob, R. A., Poudel, B., Mamtimin, M., and Zhang, J.
Integrated physics-based modeling reveals a thermodynamic gap in small modular reactor load following.
Manuscript, 2026.
```

A formal citation and DOI will be added when available.

---

## License

Please see the `LICENSE` file for terms of use.

---

## Acknowledgments

This work was supported by the U.S. Department of Energy under Award DE-NE0009296.

---

## Contact

For questions about the model or repository, please contact the repository maintainers.
