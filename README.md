# Repository for HNL analysis with highly displaced dileptons

To be cloned within an existing CMSSW installation. Tested on the base release CMSSW 13_0_13.

Contains the following projects:
* NanoProd: Production of a custom NanoAOD, including fitted secondary vertices for the dimuons.

## How to run miniAOD->nanoAOD skims production with HNLProd

### How to install
```sh
git clone --recursive git@github.com:cms-hnl/HNLProd.git
```

### Loading environment
Following command activates the framework environment:
```sh
cd HNLProd
source $PWD/env.sh nano_prod
```
### How to run miniAOD->nanoAOD skims production

Production should be run on the server that have the crab stageout area mounted to the file system.
1. Load environment
   ```sh
   source $PWD/env.sh nano_prod
   voms-proxy-init -voms cms -rfc -valid 192:00
   ```

1. Modify output and other site-specific settings in `NanoProd/crab/overseer_cfg.yaml`. In particular:
   - params/outputs/crabOutput
   - params/outputs/finalOutput

1. Test that the code works locally (take one of the miniAOD files as an input). E.g.
   ```sh
   mkdir -p tmp && cd tmp
   cmsEnv python3 $ANALYSIS_PATH/RunKit/nanoProdWrapper.py customise=HNL/NanoProd/DiMuon_cff.nanoAOD_customizeDisplacedDiMuon_Run3 maxEvents=100 sampleType=mc era=Run3_2022 inputFiles=file:/eos/user/a/amascell/HNL/output/Run3/Run3Summer22/HeavyNeutrino_trilepton_M-1.5_V-0.032_mu_massiveAndCKM_LO/miniAOD_1.root writePSet=True output="nano.root;output;$ANALYSIS_PATH/HNL/NanoProd/config/skim_mu.yaml;skim;skim_failed"
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
   python RunKit/crabOverseer.py --cfg HNL/NanoProd/crab/overseer_cfg.yaml HNL/NanoProd/crab/Run3_2022/FILE1.yaml HNL/NanoProd/crab/Run3_2022/FILE2.yaml ...
   ```
   - Use `HNL/NanoProd/crab/Run3_2022/*.yaml` to submit all the tasks
   - For more information about available command line arguments run `python RunKit/crabOverseer.py --help`
   - For consecutive runs, if there are no modifications in the configs, it is enough to run `crabOverseer` without any arguments:
     ```sh
     python RunKit/crabOverseer.py
     ```