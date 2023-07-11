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

## How to run miniAOD->nanoAOD skims production with HNLProd

### How to install
Connect to CentOS8 machine and clone the central repository.
```sh
git clone --recursive git@github.com:cms-hnl/HNLProd.git
```

### Loading environment
Following command activates the framework environment:
```sh
cd HNLProd
source $PWD/env.sh nano_prod
```
## How to run miniAOD->nanoAOD skims production

Production should be run on the server that have the crab stageout area mounted to the file system.
1. Load environment on CentOS8 machine
   ```sh
   source $PWD/env.sh nano_prod
   voms-proxy-init -voms cms -rfc -valid 192:00
   ```

1. Modify output and other site-specific settings in `config/overseer_cfg.yaml`. In particular:
   - site
   - crabOutput
   - localCrabOutput
   - finalOutput
   - renewKerberosTicket

1. Test that the code works locally (take one of the miniAOD files as an input). E.g.
   ```sh
   mkdir -p tmp && cd tmp
   cmsEnv python3 $ANALYSIS_PATH/RunKit/nanoProdWrapper.py customise=HNL/NanoProd/DiMuon_cff.nanoAOD_customizeDisplacedDiMuon_Run2 maxEvents=2000 sampleType=data storeFailed=True era=Run2_2018 inputFiles=file:/eos/cms/store/group/phys_tau/kandroso/miniAOD_UL18/SingleMuon.root skimCfg=$ANALYSIS_PATH/HNL/NanoProd/config/skim_mu.yaml writePSet=True skimSetup=skim skimSetupFailed=skim_failed createTar=False
   cmsEnv $ANALYSIS_PATH/RunKit/crabJob.sh
   ```
   Check that output file `nano_0.root` is created correctly. After that, you can remove `tmp` directory:
   ```sh
   cd $ANALYSIS_PATH
   rm -r tmp
   ```
1. Test that post-processing task is known to law:
   ```sh
   law index
   law run ProdTask --help
   ```

1. Test a dryrun crab submission
   ```sh
   python RunKit/crabOverseer.py --work-area crab_test --cfg HNL/NanoProd/crab/overseer_cfg.yaml --no-loop HNL/NanoProd/crab/crab_test.yaml
   ```
   - NB. Crab estimates of processing time will not be accurate, ignore them.
   - After the test, remove `crab_test` directory:
     ```sh
     rm -r crab_test
     ```

1. Submit tasks using `RunKit/crabOverseer.py` and monitor the process.
   It is recommended to run `crabOverseer` in screen.
   ```sh
   python RunKit/crabOverseer.py --cfg HNL/NanoProd/crab/overseer_cfg.yaml HNL/NanoProd/crab/Run2_2018/FILE1.yaml HNL/NanoProd/crab/Run2_2018/FILE2.yaml ...
   ```
   - Use `HNL/NanoProd/crab/Run2_2018/*.yaml` to submit all the tasks
   - For more information about available command line arguments run `python RunKit/crabOverseer.py --help`
   - For consecutive runs, if there are no modifications in the configs, it is enough to run `crabOverseer` without any arguments:
     ```sh
     python RunKit/crabOverseer.py
     ```
## Example of miniAOD->nanoAOD skims production for Run 3
```sh
mkdir -p tmp && cd tmp
cmsEnv python3 $ANALYSIS_PATH/RunKit/nanoProdWrapper.py maxEvents=2000 sampleType=mc era=Run3_2022 inputFiles=file:/eos/home-k/kandroso/cms-hnl/Run3/Run3Summer22/HeavyNeutrino_trilepton_M-4_V-0.00151_mu_massiveAndCKM_LO/miniAOD_1.root writePSet=True createTar=False customise=HNL/NanoProd/DiMuon_cff.nanoAOD_customizeDisplacedDiMuon_Run3 skimCfg=$ANALYSIS_PATH/HNL/NanoProd/config/skim_mu.yaml skimSetup=skim skimSetupFailed=skim_failed
cmsEnv $ANALYSIS_PATH/RunKit/crabJob.sh
```