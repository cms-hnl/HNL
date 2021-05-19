# Configurations dependent on the sample type.

import sys
from sets import Set
import FWCore.ParameterSet.Config as cms
import os

mcSampleTypes = Set([ 'MC_16', 'MC_17', 'MC_18' ])
dataSampleTypes = Set([ 'Run2016',  'Run2016UL', 'Run2017', 'Run2018ABC', 'Run2018D' ])

periodDict = { 'MC_16' : 'Run2016',
               'MC_16_UL' : 'Run2016',
               'Run2016' : 'Run2016',
               'Run2016UL' : 'Run2016',
               'MC_17' : 'Run2017',
               'Run2017' : 'Run2017',
               'MC_18' : 'Run2018',
               'Run2018ABC' : 'Run2018',
               'Run2018D' : 'Run2018',
             }

globalTagMap = { 'MC_16' : '102X_mcRun2_asymptotic_v7',
                 'Run2016' : '102X_dataRun2_v12',
                 'Run2016UL' : '106X_dataRun2_v32',
                 'MC_17' : '102X_mc2017_realistic_v7',
                 'Run2017' : '102X_dataRun2_v12',
                 'MC_18' : '102X_upgrade2018_realistic_v21',
                 'Run2018ABC' : '102X_dataRun2_v12',
                 'Run2018D' : '102X_dataRun2_Prompt_v15',
               }

def isData(sampleType):
    isData = sampleType in dataSampleTypes
    if not isData and not sampleType in mcSampleTypes:
        raise RuntimeError("ERROR: unknown sample type = '{}'".format(sampleType))
    return isData

def getPeriod(sampleType):
    if sampleType not in periodDict:
        raise RuntimeError("ERROR: unknown sample type = '{}'".format(sampleType))
    return periodDict[sampleType]

def getGlobalTag(sampleType):
    if sampleType not in globalTagMap:
        raise RuntimeError("ERROR: unknown sample type = '{}'".format(sampleType))
    return globalTagMap[sampleType]

def createProcess(processName, sampleType):
    period = getPeriod(sampleType)
    if period == 'Run2018':
        from Configuration.StandardSequences.Eras import eras
        return cms.Process(processName, eras.Run2_2018, eras.run2_nanoAOD_102Xv1)
    elif period == 'Run2016':
        from Configuration.StandardSequences.Eras import eras
        return cms.Process(processName, eras.Run2_2016, eras.run2_nanoAOD_106Xv1)
    else:
        raise RuntimeError('Period = "{}" is not supported.'.format(period))
