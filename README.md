# Repository for HNL analysis with highly displaced dileptons

To be cloned within an existing CMSSW installation. Tested on the 10_2_X branch (base release CMSSW 10_2_22).

Contains the following projects:
* NanoProd: Production of a custom NanoAOD, including fitted secondary vertices for the dimuons.

## How to run a private MC production

- currently using centrally produced grid packs used for the official production
- reference scripts to run 3 production steps (gen-sim, premix-hnl and reco): [run_hnl_prod.sh](https://github.com/cms-hnl/HNL/blob/master/RecoProd/scripts/run_hnl_prod.sh) and [run_hnl_prod_loop.sh](https://github.com/cms-hnl/HNL/blob/master/RecoProd/scripts/run_hnl_prod_loop.sh)

### Installation
```shell
cmsrel CMSSW_10_2_25
cd CMSSW_10_2_25/src
cmsenv
git clone git@github.com:cms-hnl/HNL.git
```

### How to run run_hnl_prod.sh
- Before running, variables defined in the beginning of run_hnl_prod.sh should be set accordingly:
  - __FLAV__ - flavour (e, mu, tau)
  - __MASS__ - HNL mass
  - __OUT__ - directory where output root files will be placed
- Make sure that a valid grid certificate is initialised (needed to access `Neutrino_E-10_gun` dataset for PU mixing)
  ```shell
  voms-proxy-init -rfc -voms cms -valid 192:00
  ```

```shell
./HNL/RecoProd/scripts.run_hnl_prod.sh SEED EVENTS
```
where
- __SEED__ - a positive integer number to initiate the random generator for event generation
- __EVENTS__ - number of events to generate. Current implementation works only for EVENTS <= 5000

The final output file will be located in `$OUT/HeavyNeutrino_trilepton_M-${MASS}_${FLAV}/HeavyNeutrino_trilepton_M-${MASS}_${FLAV}_S${SEED}.root`

`run_hnl_prod_loop.sh` can be used to consecutively run production for a range of seeds.
