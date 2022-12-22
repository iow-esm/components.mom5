# Purpose, Description

This is the Modular ocean model (MOM5) component of the IOW ESM.
The starting point for the IOW ESM development was a IOW-adapted version of the code that is originating from
the following commit on https://github.com/mom-ocean/MOM5

Autor:			Nicholas Hannah <nicholash@users.noreply.github.com>
Datum:			vor 7 Jahren (24.03.2014 04:17:14)
Commit Hash:	f406b4c5b4bbece3b0ae7f376a4ba90ea68ffb1e
Kind:			Commit Index

Merge pull request #33 from nicholash/mom-87-various-bugs

Mom 87 various bugs

Please refer to https://mom-ocean.github.io/docs/userguide/ for more information.


# Authors

* SK      (sven.karsten@io-warnemuende.de)
* HR      (hagen.radtke@io-warnemuende.de)


# Versions

## 1.01.00 (latest release)

| date        | author(s)   | link      |
|---          |---          |---        |
| 2022-12-22  | SK          | [1.01.00](https://git.io-warnemuende.de/iow_esm/components.mom5/src/branch/1.01.00)  |  

<details> 

### changes
* updated ERGOM module to newest available version
  * Caution: this code is not compatible with current 8nm example setups!
  * 8nm examples have to be updated as well
* build with compiler flags that lead to more stablility (not so many NaN problems)

### dependencies
* OASIS3-MCT libraries
* see build scripts for more dependencies  

### known issues
* in coupled model this version leads to too cold winter temperatures 
  when evaluated from 1979-2009
* model is not running on phy-2

### tested with
* intensively tested on both HLRN machines
  * using example setups available under:
    (coupled) /scratch/usr/mviowmod/IOW_ESM/setups/
              MOM5_Baltic-CCLM_Eurocordex/example_3nm_0.22deg/1.00.00
    (uncoupled) /scratch/usr/mviowmod/IOW_ESM/setups/
                MOM5_Baltic/example_3nm/1.00.00
   * results exhibit known issues
  * can be built and run on Haumea but output is not intensively tested


## 1.00.00 

| date        | author(s)   | link      |
|---          |---          |---        |
| 2022-01-28  | SK, HR      | [1.00.00](https://git.io-warnemuende.de/iow_esm/components.mom5/src/branch/1.00.00)  |  

<details> 

### changes
* initial release
* main codes is based on branch 5.1.0 on https://github.com/mom-ocean/MOM5.git 
* coupling interface has been implemented by authors
* adapted Fopts and added build scripts for 
  * both HLRN supercomputers
  * Rostock University's cluster Haumea
  * for IOW phy-2 machine (model is not running here)
* added OASIS interface to flux calculator
  * send SST, fraction of ice and albedo
  * send variables to calculate fluxes for evaporation, latent and sensible heat,
    momentum 
  * receive radiation and precipitation
  * receive calculated fluxes in corresponding variables that have been fed by forcing data in the uncoupled mode
* executable can run in coupled model or as stand alone
  * use parameter `type_atmos='flux_calculator'` in `coupler_nml` in `input.nml` for control, see also documentation in
      iow_esm/main project

### dependencies
* OASIS3-MCT libraries
* see build scripts for more dependencies  

### known issues
* in coupled model this version leads to too cold winter temperatures 
  when evaluated from 1979-2009
* model is not running on phy-2

### tested with
* intensively tested on both HLRN machines
  * using example setups available under:
    (coupled) /scratch/usr/mviowmod/IOW_ESM/setups/
              MOM5_Baltic-CCLM_Eurocordex/example_8nm_0.22deg/1.00.00
    (uncoupled) /scratch/usr/mviowmod/IOW_ESM/setups/
                MOM5_Batlic/example_8nm/1.00.00
   * results exhibit known issues
  * can be built and run on Haumea but output is not intensively tested