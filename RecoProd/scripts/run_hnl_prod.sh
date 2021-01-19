#!/bin/bash

FLAV=mu
MASS=4
#FRAGMENT=Configuration/GenProduction/python/EXO-RunIIFall18wmLHEGS-00886-fragment.py
COND=102X_upgrade2018_realistic_v15
BS=Realistic25ns13TeVEarly2018Collision
GEOM=DB:Extended
ERA=Run2_2018
EVENTS=$2
N_THREADS=6
SEED=$1
OUT=output
CFG=cfg
PU="dbs:/Neutrino_E-10_gun/RunIISummer17PrePremix-PUAutumn18_102X_upgrade2018_realistic_v15-v1/GEN-SIM-DIGI-RAW"

NAME=HeavyNeutrino_trilepton_M-${MASS}_${FLAV}_S${SEED}

FRAGMENT=Configuration/GenProduction/python/HeavyNeutrino_trilepton_M-${MASS}_${FLAV}-fragment.py
if ! [ -f $FRAGMENT ]; then
    echo "Fragment $FRAGMENT not found. Creating it..."
    GRIDPACK="$(ls /cvmfs/cms.cern.ch/phys_generator/gridpacks/slc6_amd64_gcc481/13TeV/madgraph/V5_2.4.2/exo_heavyNeutrino/displaced_2017/v1/HeavyNeutrino_trilepton_M-${MASS}_V-*_${FLAV}_massiveAndCKM_LO_slc6_amd64_gcc481_CMSSW_7_1_30_tarball.tar.xz)"
    N_GRIDPACK=$(echo "$GRIDPACK" | wc -l)
    if [ $N_GRIDPACK -ne 1 ]; then
        echo "Multiple gridpacks were found."
        exit 1
    fi
    if ! [ -f $GRIDPACK ]; then
        echo "Grid pack not found."
        exit 1
    fi
    mkdir -p Configuration/GenProduction/python/
    cat <<EOF > $FRAGMENT
import FWCore.ParameterSet.Config as cms

externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    args = cms.vstring('$GRIDPACK'),
    nEvents = cms.untracked.uint32(5000),
    numberOfParameters = cms.uint32(1),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')
)

# Link to cards:
# https://github.com/cms-sw/genproductions/tree/e3bf0a9b8180b78938a34e4d06036d5dd096d8ef/bin/MadGraph5_aMCatNLO/cards/production/2017/13TeV/exo_heavyNeutrino_LO/

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *
from Configuration.Generator.PSweightsPythia.PythiaPSweightsSettings_cfi import *

generator = cms.EDFilter("Pythia8HadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        pythia8PSweightsSettingsBlock,
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CP5Settings',
                                    'pythia8PSweightsSettings',
                                    )
    )
)
EOF
    scram b
fi

CURRENT_OUTPUT=None

function run_cmd {
    echo "> $@"
    "$@"
    RESULT=$?
    if [ $RESULT -ne 0 ] ; then
        echo "Error while rinning '$@'"
        if [ -f $CURRENT_OUTPUT ] ; then
            echo "Removing current output $CURRENT_OUTPUT"
            rm "$CURRENT_OUTPUT"
        fi
        exit 1
    fi
}

if [ -f $OUT/HeavyNeutrino_trilepton_M-${MASS}_${FLAV}/$NAME.root ]; then
    echo "$OUT/HeavyNeutrino_trilepton_M-${MASS}_${FLAV}/$NAME.root already exists."
    exit 0
fi

if [ -f $OUT/interm_steps/${NAME}_SIM.root ]; then
    echo "$OUT/interm_steps/${NAME}_SIM.root already exists. Continuing to the next step."
else
    CURRENT_OUTPUT=$OUT//interm_steps/${NAME}_SIM.root
    run_cmd cmsDriver.py $FRAGMENT --python_filename $CFG/${NAME}_SIM_cfg.py --eventcontent RAWSIM,LHE \
        --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM,LHE \
        --fileout file:$OUT/interm_steps/${NAME}_SIM.root --conditions $COND --beamspot $BS \
        --customise_commands process.RandomNumberGeneratorService.externalLHEProducer.initialSeed="int(${SEED})" \
        --step LHE,GEN,SIM --geometry $GEOM --era $ERA --no_exec --mc --nThreads $N_THREADS -n $EVENTS

    run_cmd cmsRun $CFG/${NAME}_SIM_cfg.py
fi

if [ -f $OUT/interm_steps/${NAME}_PREMIX.root ]; then
    echo "$OUT/interm_steps/${NAME}_PREMIX.root already exists. Continuing to the next step."
else
    CURRENT_OUTPUT=$OUT//interm_steps/${NAME}_PREMIX.root
    run_cmd cmsDriver.py --python_filename $CFG/${NAME}_PREMIX_cfg.py --eventcontent PREMIXRAW \
        --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-RAW \
        --fileout file:$OUT/interm_steps/${NAME}_PREMIX.root --pileup_input $PU --conditions $COND \
        --step DIGI,DATAMIX,L1,DIGI2RAW,HLT:@relval2018 --procModifiers premix_stage2 --geometry $GEOM \
        --filein file:$OUT/interm_steps/${NAME}_SIM.root --datamix PreMix --era $ERA --no_exec --mc \
        --nThreads $N_THREADS -n $EVENTS

    run_cmd cmsRun $CFG/${NAME}_PREMIX_cfg.py
fi

CURRENT_OUTPUT=$OUT/HeavyNeutrino_trilepton_M-${MASS}_${FLAV}/${NAME}.root
run_cmd cmsDriver.py --python_filename $CFG/${NAME}_RECO_cfg.py --eventcontent FEVTDEBUGHLT \
    --customise Configuration/DataProcessing/Utils.addMonitoring --datatier RECOSIM \
    --fileout file:$OUT/HeavyNeutrino_trilepton_M-${MASS}_${FLAV}/$NAME.root --conditions $COND \
    --step RAW2DIGI,L1Reco,RECO,RECOSIM,EI --procModifiers premix_stage2 \
    --filein file:$OUT/interm_steps/${NAME}_PREMIX.root --era $ERA --runUnscheduled \
    --no_exec --mc --nThreads $N_THREADS -n $EVENTS

run_cmd cmsRun $CFG/${NAME}_RECO_cfg.py
