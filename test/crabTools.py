#!/usr/bin/env python
import os
import sys
import ROOT
import subprocess
from optparse import OptionParser
import optparse

# General campaign parameters
PWD = os.getenv('PWD') + '/'
campaign = 'UL18'
outputDir = '/store/group/phys_btag/performance/' + campaign + '/'
sampleListDir = './samples/' + campaign + '/'
maxFilesPerRange = 600
splitFilesPerRange = 500
treeHeaderName = 'RunIISummer19UL18MiniAOD-106X'
storageSite = 'T2_CH_CERN'
storageDir = '/eos/cms'
storageAddress = 'root://eoscms.cern.ch/' + storageDir
#lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
cfgDefaults = '2018_UltraLegacy'
datasetDict = { 'ttbar'  : { 'datatype'   : 'MC',
                             'directory'  : 'TTToHadronic_TuneCP5_13TeV-powheg-pythia8_RunIISummer19UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1',
                             'parameters' : '',
                             'dataname'   : '',
                             'rangename'  : '',
                             'xsecname'   : '' },
                'QCD'    : { 'datatype'   : 'MC',
                             'directory'  : 'QCD_TuneCP5_13TeV_pythia8_RunIISummer19UL17MiniAOD-106X_mc2017_realistic_v6',
                             'parameters' : ', \'producePtRelTemplate=True\'',
                             'dataname'   : 'MCInclusive', 
                             'rangename'  : 'PtHat',
                             'xsecname'   : 'Inclusive' }, 
                'QCDMu'  : { 'datatype'   : 'MC', 
                             'directory'  : 'QCD_MuEnrichedPt5_TuneCP5_13TeV_pythia8_RunIISummer19UL17MiniAOD-106X_mc2017_realistic_v6',
                             'parameters' : '',
                             'dataname'   : 'MonteCarlo', 
                             'rangename'  : 'PtHat',
                             'xsecname'   : '' },
                'JetHT'  : { 'datatype'   : 'Data', 
                             'directory'  : 'JetHT_Run2017-09Aug2019_UL2017',
                             'parameters' : ', \'producePtRelTemplate=True\'',
                             'dataname'   : 'Jet', 
                             'rangename'  : 'Run' },
                'BTagMu' : { 'datatype'   : 'Data', 
                             'directory'  : 'BTagMu_Run2017-09Aug2019_UL2017',
                             'parameters' : '',
                             'dataname'   : 'BTagMu', 
                             'rangename'  : '' },
}

# Method to write crab cfg file
def writeConfigurationFile(crabDir, dataset, sample, dasDataset):

    fCfg = open('./' + crabDir + '/BTA_cfg.py','w')

    fCfg.write('from WMCore.Configuration import Configuration\n')

    fCfg.write('config = Configuration()\n')

    fCfg.write('config.section_(\'General\')\n')
    fCfg.write('config.General.requestName = \'' + sample + '\'\n')
    #fCfg.write('config.General.workArea = \'' + crabDir + '\'\n')
    fCfg.write('config.General.transferOutputs = True\n')
    fCfg.write('config.General.transferLogs = False\n')

    fCfg.write('config.section_(\'JobType\')\n')
    fCfg.write('config.JobType.pluginName = \'Analysis\'\n')    
    fCfg.write('config.JobType.psetName = \'' + PWD + 'runBTagAnalyzer_cfg.py\'\n')
    cfgParams = '\'runOnData=False\', \'defaults=' + cfgDefaults + '\'' + datasetDict[dataset]['parameters'] + ', \'groups=Luca\''
    if 'defaults' in samples[sample]:
        cfgParams = cfgParams.replace(cfgDefaults, cfgDefaults+samples[sample]['defaults'])
    if datasetDict[dataset]['datatype']=='Data':
        cfgParams = cfgParams.replace('runOnData=False', 'runOnData=True')
    fCfg.write('config.JobType.pyCfgParams = [' + cfgParams + ']\n')
    #fCfg.write('config.JobType.inputFiles = [\'\']\n')
    #fCfg.write('config.JobType.outputFiles = [\'\']\n')
    fCfg.write('config.JobType.allowUndistributedCMSSW = True\n')
    fCfg.write('config.JobType.sendExternalFolder = True\n')

    fCfg.write('config.section_(\'Data\')\n')
    fCfg.write('config.Data.inputDataset = \'' + dasDataset + '\'\n')
    fCfg.write('config.Data.allowNonValidInputDataset = True\n')
    fCfg.write('config.Data.inputDBS = \'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/\'\n')
    if datasetDict[dataset]['datatype']=='MC':
        fCfg.write('config.Data.splitting = \'FileBased\'\n')
        fCfg.write('config.Data.unitsPerJob = 1\n')
    else:
        fCfg.write('config.Data.splitting = \'LumiBased\' # Recommended for data\n') 
        fCfg.write('config.Data.unitsPerJob = 100\n')
        fCfg.write('config.Data.lumiMask = \'' + lumiMask + '\'\n')
    fCfg.write('config.Data.outLFNDirBase = \'' + outputDir + datasetDict[dataset]['datatype'] + '/\'\n')
    fCfg.write('config.Data.publication = False\n')
    fCfg.write('config.Data.ignoreLocality = False\n')
               
    fCfg.write('config.section_(\'User\')\n')

    fCfg.write('config.section_(\'Site\')\n')
    fCfg.write('config.Site.storageSite = \'' + storageSite + '\'\n')
    #fCfg.write('config.Site.blacklist = [ ]\n')
    #fCfg.write('config.Site.whitelist = [ ]\n')

    fCfg.close() 

# Method to write tree header file
def writeTreeHeaderFile(datasets):

    cleanCrab = True

    fTree = open('./' + sampleListDir + treeHeaderName + '.h','w')

    fTree.write('TString TreeContentFlag = "";\n\n')

    for dataset in datasets:

        isMC = datasetDict[dataset]['datatype']=='MC'

        fTree.write('// ' + dataset + ' ' + datasetDict[dataset]['datatype'] + '\n')
        fTree.write('TString EOSPath' + dataset + ' = "' + storageAddress + outputDir + datasetDict[dataset]['datatype'] + '/' + datasetDict[dataset]['directory'] + '/";\n')

        samples = {}
        if os.path.exists(sampleListDir + dataset + '.py'):
            exec(open(sampleListDir + dataset + '.py').read())
        else:
            print 'Sample list for', dataset, 'dataset not found -> crab directories will not be deleted' 
            cleanCrab = False
            fTree.write('\n')
            continue

        nSamples = len(samples)

        sampleRanges = '{'
        samplenTrees = '{'
        sampleEvents = '{'
        sampleCrSecs = '{'
        for sample in samples:
            sampleRange = sample.split('_', 1)[1].replace('Pt_', 'Pt-')
            nTrees = '0'
            events = '0'
            if os.path.exists('crab_' + dataset + '/done_' + sample):
                outputSampleDir = storageDir + outputDir + datasetDict[dataset]['datatype'] + '/' + datasetDict[dataset]['directory'] + '/' + sampleRange
                nTrees = subprocess.check_output('ls ' + outputSampleDir + '/*.root | grep -c root', shell=True).strip('\n')
                if isMC:
                    chain = ROOT.TChain('btagana/ttree')
                    chain.Add(outputSampleDir + '/*.root')
                    print chain.GetEntries()
                    events = str(chain.GetEntries())
            else:
                print 'Tree production for', sample, 'sample not completed -> crab directories will not be deleted' 
                cleanCrab = False
            nSplitSample = 1 
            if int(nTrees)>maxFilesPerRange:
                nSplitSample = (int(nTrees)-1)/splitFilesPerRange + 1
                nSamples += nSplitSample - 1
            for sampleSplit in range(1, nSplitSample+1):
                splitPart = '' if nSplitSample==1 else ':'+str(sampleSplit)+'of'+str(nSplitSample)
                sampleRanges += ' "' + sampleRange + splitPart + '",'
                samplenTrees += ' ' + nTrees + ','
                if isMC:
                    sampleCrSecs += ' ' + samples[sample]['crossSection'] + ','
                    sampleEvents += ' ' + events + ','

        datasetRangeFlag = datasetDict[dataset]['dataname'] + datasetDict[dataset]['rangename'] + 'Range'
        fTree.write('const int n' + datasetRangeFlag + 's = ' + str(nSamples) + ';\n')

        Name = '' if isMC else 'Name'
        fTree.write('TString ' + datasetRangeFlag + Name + '[n' + datasetRangeFlag + 's] = ' + sampleRanges.strip(',') + ' };\n')
        fTree.write('int n' + datasetDict[dataset]['dataname'] + 'Trees[n' + datasetRangeFlag + 's] = ' + samplenTrees.strip(',') + ' };\n')

        if isMC:
            fTree.write('double CrossSection' + datasetDict[dataset]['xsecname'] + '[n' + datasetRangeFlag + 's] = ' + sampleCrSecs.strip(',') + ' };\n')
            fTree.write('double GeneratedEvents' + datasetDict[dataset]['xsecname'] + '[n' + datasetRangeFlag + 's] = ' + sampleEvents.strip(',') + ' };\n')

        fTree.write('\n')

    fTree.write('float PileUpScenario[] = { }\n')

    fTree.close()

    return cleanCrab

# Main
if __name__ == '__main__':

    # Input parameters
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)

    parser.add_option('-d', '--dataset', dest='dataset', help='Dataset to be processed' , default='all')
    parser.add_option('-T', '--seldata', dest='seldata', help='Select samples'          , default='all')
    parser.add_option('-E', '--excdata', dest='excdata', help='Exclude samples'         , default='none')
    parser.add_option('-a', '--action',  dest='action',  help='Action to perform'       , default='status')
    (opt, args) = parser.parse_args()

    datasets = datasetDict.keys() if opt.dataset=='all' else opt.dataset.split(',') 

    if opt.action=='terminate':

        if writeTreeHeaderFile(datasets):
            for dataset in datasetDict:
                #os.system('rm -rf crab_' + dataset)
                print 'Crab directory for', dataset, 'dataset deleted'

        else:
            print 'Skipping crab directory deletion'

        exit()

    seldata = opt.seldata.split(',')
    excdata = opt.excdata.split(',')

    actions = [ ]

    if opt.action=='update': opt.action = 'status,resubmit,store'
    if opt.action=='clean': opt.action = 'kill,delete'

    for action in [ 'submit', 'status', 'resubmit', 'kill', 'store', 'delete' ]:
 	if action in opt.action.split(','):
	    actions.append(action)

    for action in opt.action.split(','):
        if action not in actions:
	    print 'Unknown action:', action, ' --> exiting'
            exit()

    for dataset in datasets:

        print ''

        if dataset not in datasetDict:
            print dataset + ': dataset not supported --> ignoring'
            continue

        crabDir = 'crab_' + dataset
        if 'submit' in actions:
            os.system('mkdir -p ' + crabDir)
        elif not os.path.isdir(crabDir):
            print dataset + ': crab directory ' + crabDir + ' not existing --> skipping'
            continue

        outputCrabDir = storageDir + outputDir + datasetDict[dataset]['datatype'] + '/'
        if not os.path.isdir(outputCrabDir):
            if 'submit' in actions:
                os.system('mkdir -p ' + outputCrabDir)
            else:
                print dataset + ': output crab directory ' + outputCrabDir + ' not existing --> skipping'
                continue
                
        outputDatasetDir = outputCrabDir + datasetDict[dataset]['directory'] + '/'
        os.system('mkdir -p ' + outputDatasetDir)

        print dataset + ': looping on samples'

        samples = {}
        if os.path.exists(sampleListDir + dataset + '.py'):
            exec(open(sampleListDir + dataset + '.py').read())
        else:
            print 'Sample list for', dataset, 'dataset not found --> skipping'
            continue

        lengthSampleName = max(len(sample) for sample in samples)

        for sample in samples:
            if (opt.seldata=='all' or sample in seldata) and sample not in excdata:

                crabSampleDir = crabDir + '/crab_' + sample
                outputSampleDir = outputDatasetDir + sample.split('_', 1)[1].replace('Pt_', 'Pt-') + '/'

                taskName = ''
                hasFailed = True
                nFinished = 1

                header = [ '  Sample ' + sample + ':' + ' '*(lengthSampleName-len(sample)) ]
                header.append(' '*len(header[0])) 

                if os.path.exists(crabSampleDir.replace('/crab_', '/done_')):
                    if 'delete' in actions:
                        os.system('rm -r ' + outputSampleDir)
                        os.system('rm -r ' + crabSampleDir.replace('/crab_', '/done_'))
                        print header[0], 'deleted'
                    else:
                        ntrees = int(subprocess.check_output('ls ' + outputSampleDir + '*.root | grep -c root', shell=True))
                        print header[0], ntrees, 'trees done'
                    continue

                for action in actions:
                    
                    if action=='delete':
                        pass

                    elif action=='submit':	
                
                        if os.path.exists(crabSampleDir + '/.requestcache'):
                            print header[0], 'crab job already submitted'

                        else:

                            if os.path.isdir(crabSampleDir):
                                os.system('rm -r ' + crabSampleDir)

                            print header[0], 'submitting crab job'
                            writeConfigurationFile(crabDir, opt.dataset, sample, samples[sample]['das'])
                            os.system('cd ' + crabDir + '; crab submit BTA_cfg.py; rm BTA_cfg.py; rm BTA_cfg.pyc; cd -')
                            #submit_result = subprocess.check_output('cd ' + crabDir + '; crab submit BTA_cfg.py; rm BTA_cfg.py; rm BTA_cfg.pyc; cd -', shell=True)
                            print header[1], 'crab job submitted'

                    else:

                        if not os.path.isdir(crabSampleDir):
                            print header[0], 'crab job not submitted'
                            break

                        elif not os.path.exists(crabSampleDir + '/.requestcache'):
                            print header[0], 'cannot find .requestcache file in CRAB project'   
                            break

                        else:
                        
                            if action=='resubmit' and not hasFailed:
                                print header['status' in actions], 'there are not failed jobs --> will not resubmit'
                                continue

                            if action=='store' and nFinished==0: 
                                print header['status' in actions], 'not all jobs done --> will not store'
                                break

                            crabAction = action if action!='store' else 'status'

                            try:
                                
                                if action!='store' or 'status' not in actions:
                                    action_result = subprocess.check_output('crab ' + crabAction + ' -d ' + crabSampleDir, shell=True)

                                if action=='status':

                                    hasFailed = False
                                    nFinished = 0

                                    idx = 0 
                                    for job_status in [ 'finished', 'transferring', 'running' , 'toRetry', 'failed', 'idle', 'unsubmitted' ]:
                                        for line in action_result.split('\n'):
                                            if job_status in line and '%' in line and '(' in line and '/' in line and ')' in line and 'jobs' not in line:

                                                print header[idx>0], job_status, line.split(job_status)[1]
                                                if job_status=='finished' and '100.0%' in line:
                                                    nFinished = int(line.split('(')[1].split('/')[0])
                                                if job_status=='failed': hasFailed = True
                                                idx += 1

                                    taskName = [line for line in action_result.split('\n') if 'Task name:' in line][0].split('\t')[3]

                                else:
 
                                    if action=='store':

                                        isStored = False

                                        if 'status' in actions:
                                            pass

                                        else:

                                            nFinished = 0
                                            for line in action_result.split('\n'):
                                                if 'finished' in line and '100.0%' in line:
                                                    nFinished = int(line.split('(')[1].split('/')[0])
                                            if nFinished==0:
                                                print header['status' in actions], 'not all jobs done --> will not store'
                                                break
                                            taskName = [line for line in action_result.split('\n') if 'Task name:' in line][0].split('\t')[3]
                                
                                        outputJobDir = outputCrabDir + samples[sample]['das'].split('/')[1] + '/crab_' + sample + '/' + taskName.split(':')[0] + '/000*/'

                                        if subprocess.check_output('ls ' + outputJobDir, shell=True)=='':
                                            print header['status' in actions], 'output job directory ' + outputJobDir + ' empty --> will not store'
                                            break

                                        if nFinished!=int(subprocess.check_output('ls ' + outputJobDir + '*.root | grep -c root', shell=True)):
                                            print header['status' in actions], 'number of trees does not match number of jobs --> will not store'
                                            break
                                            
                                        os.system('mkdir -p ' + outputSampleDir)
                                        
                                        if subprocess.check_output('ls ' + outputSampleDir, shell=True)!='':
                                            print header['status' in actions], 'output directory ' + outputSampleDir + ' not empty --> will not store'
                                            break

                                        os.system('mv ' + outputJobDir + '/*.root ' + outputSampleDir)
                                        
                                        if subprocess.check_output('ls ' + outputSampleDir, shell=True)=='':
                                            print header['status' in actions], 'output directory ' + outputSampleDir + ' still empty --> please check'
                                            break

                                        if nFinished!=int(subprocess.check_output('ls ' + outputSampleDir + '*.root | grep -c root', shell=True)):
                                            print header['status' in actions], 'number of stored trees does not match number of jobs --> please check'
                                            break

                                    if action=='store' or action=='kill':
                                            
                                        subSampleDir = '/' if datasetDict[dataset]['datatype']=='MC' else '/crab_' + sample
                                        os.system('rm -r ' + crabSampleDir)
                                        os.system('rm -r ' + outputCrabDir + samples[sample]['das'].split('/')[1] + subSampleDir)

                                        if action=='store': 
                                            os.system('touch ' + crabSampleDir.replace('/crab_', '/done_'))

                                        elif action=='kill':
                                            os.system('rm -r -f ' + outputSampleDir)
                                            
                                    print header['status' in actions], action.replace('mit', 'mitt').replace('ore', 'or') + 'ed'


                            except subprocess.CalledProcessError as e:
                                raise Warning("command '{}' return with errorclean (code {}): {}".format(e.cmd, e.returncode, e.output))          
	
                    



            
    
