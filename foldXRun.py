#!/usr/bin/env python
from fabric.api import local, env, settings, lcd, hide
import re
from tempfile import mkdtemp
from shutil import rmtree, copy
from os.path import isfile, join, normpath, exists, abspath, dirname, basename
import sys
from os import listdir, makedirs
import numpy
import csv

defaultOptions = """\
<OPTIONS>FOLDX_optionfile;
<Temperature>298;
<R>#;
<pH>7;
<IonStrength>0.050;
<water>-CRYSTAL;
<metal>-CRYSTAL;
<VdWDesign>2;
<OutPDB>true;
<pdb_hydrogens>false;
<END>#;"""
defaultJob = """\
<JOBSTART>#;
<PDBS>%(pdbsString)s;
<BATCH>#;
<COMMANDS>FOLDX_commandfile;
%(commandsString)s
<END>#;
%(optionsString)s
<JOBEND>#;"""
defaultRunFile = """\
<TITLE>FOLDX_runscript;
%s
<ENDFILE>#;"""

def foldXRun(seqDirectory, mutDirectory, nativePDB, foldXPath = '.', outputDir = '.'):
# First make sure data directory exists and is writeable
    dataDir = abspath(join(outputDir, 'data'))
    scratchDir = abspath(join(outputDir, 'scratch'))
    foldXPath = abspath(foldXPath)
    if not exists(dataDir):
        try:
            makedirs(dataDir)
        except:
            exit("Couldn't create destination directories and they don't exist!")
    if not exists(scratchDir):
        try:
            makedirs(scratchDir)
        except:
            exit("Couldn't create scratch directory and it doesn't exist")


    results = {}
    results['seqDir'] = seqDirectory
    results['mutDir'] = mutDirectory
    results['dataDir'] = abspath(dataDir)
    # First we handle all model-model comparisons
    # Generates a list of all model[0-9].pdb in seqDirectory
    seqModels = [f for f in listdir(seqDirectory) if re.search("model[0-9]*.pdb", f) and isfile(join(seqDirectory, f))]
    mutModels = [f for f in listdir(mutDirectory) if re.search("model[0-9]*.pdb", f) and isfile(join(mutDirectory, f))]
    results['numSeqModels'] = len(seqModels)
    results['numMutModels'] = len(mutModels)
    results['nativePDB'] = nativePDB
    # We copy all model files to our scratch directory so we can work on em
    for modelFile in seqModels:
        copy(join(seqDirectory, modelFile), join(scratchDir, "seq_" + modelFile))
    for modelFile in mutModels:
        copy(join(mutDirectory, modelFile), join(scratchDir, "mut_" + modelFile))
    copy(nativePDB, join(scratchDir, "native.pdb"))
    copy(join(dirname(foldXPath), 'rotabase.txt'), scratchDir)

    # Now we setup the repair job
    repairDict = {'optionsString':defaultOptions}

    # Setup repair PDBs
    repairJobPDBs = ""
    for modelFile in seqModels:
        repairJobPDBs += "seq_%s," % modelFile
    for modelFile in mutModels:
        repairJobPDBs += "mut_%s," % modelFile
    repairJobPDBs += "native.pdb"
    repairDict['pdbsString'] = repairJobPDBs

    # Setup repair commands

    repairCommand = "<RepairPDB>#;"
    repairDict['commandsString'] = repairCommand

    # Finalize repair job.
    repairJob = defaultJob % repairDict

    # Now we handle the stability job in the same way.

    stabilizeDict = { 'optionsString':defaultOptions }

    # Setup stabilize PDBs
    stabilizeJobPDBs = ""
    for modelFile in seqModels:
        stabilizeJobPDBs += "RepairPDB_seq_%s," % modelFile
    for modelFile in mutModels:
        stabilizeJobPDBs += "RepairPDB_mut_%s," % modelFile
    stabilizeJobPDBs += "native.pdb"
    stabilizeDict['pdbsString'] = stabilizeJobPDBs

    # Setup stabilize commands

    stabilizeCommand = "<Stability>%s;" % join(dataDir, 'energyOutput.txt')
    stabilizeDict['commandsString'] = stabilizeCommand

    # Finalize stabilize job.
    stabilizeJob = defaultJob % stabilizeDict
    allJobs = repairJob + "\n" + stabilizeJob
    repairRunFileText = defaultRunFile % repairJob
    stabilizeRunFileText = defaultRunFile % stabilizeJob
    with open(join(scratchDir, 'repairRunFile.txt'), 'w') as repair, open(join(scratchDir, 'stabilizeRunFile.txt'), 'w') as stabilize:
        repair.write(repairRunFileText)
        stabilize.write(stabilizeRunFileText)
    with lcd(scratchDir), hide('everything'):
        print "Running repair job."
        local("%s -runfile repairRunFile.txt" % foldXPath)
        print "Running stabilize job."
        local("%s -runfile stabilizeRunFile.txt" % foldXPath)
        print "Done~"

    print results
    print "FoldX Complete - Now running analysis on output"
    reportAnalyze(join(dataDir, 'energyOutput.txt'), join(dataDir, 'report.txt'), basename(normpath(seqDirectory)))

def reportAnalyze(reportPath, outputPath, reportName):
    native = None
    seqs = []
    muts = []
    seqsSTDs = {}
    mutsSTDs = {}
    seqsMeans = {}
    mutsMeans = {}
    # Reading in values from the Stabilized report
    with open(reportPath, 'r') as f:
        for i in range(8):
            f.readline()
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if re.match(".*mut", row['']):
               muts.append(row)
            elif re.match(".*seq", row['']):
                seqs.append(row)
            else:
                native = row

    # Python is really bad!!
    for i in seqs:
        i['name'] = i['']
        del i['']
    for i in muts:
        i['name'] = i['']
        del i['']
    native['name'] = native['']
    del native['']

    # Calculate means, stds for each value
    for key in seqs[0]:
        # We don't want stats on these!
        if key in ('name', 'Number of Residues'):
            continue
        seqsList = [float(d[key]) for d in seqs]
        mutsList = [float(d[key]) for d in muts]

        seqsSTDs[key] = numpy.std(seqsList)
        mutsSTDs[key] = numpy.std(mutsList)
        seqsMeans[key] = numpy.mean(seqsList)
        mutsMeans[key] = numpy.mean(mutsList)

    # We have our results, now let's write them prettily.
    with open(outputPath, 'w') as f:


        # First generate a list of lines we are going to write to the file.
        lines = []
        lines.append("""\
Analysis Output for %s\n\n""" % reportName)
        lines.append("Mutant Predicted Data:\n\n")
        for key in mutsSTDs:
            lines.append("%s Mean: %f\n" % (key, mutsMeans[key]))
            lines.append("%s STD: %f\n" % (key, mutsSTDs[key]))
        lines.append("\n\n")
        lines.append("WT Predicted Data:\n\n")
        for key in mutsSTDs:
            lines.append("%s Mean: %f\n" % (key, seqsMeans[key]))
            lines.append("%s STD: %f\n" % (key, seqsSTDs[key]))
        lines.append("\n\n")
        lines.append("WT Experimental Data:\n\n")
        # We want to print the same data, but there are obviously not means or STD deviations for our single native protein.
        for key in mutsSTDs:
            lines.append("%s: %s\n" % (key, native[key]))

        # Now write!
        f.writelines(lines)



if __name__ == '__main__':
    try:
        foldXRun(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    except:
        print """\
Incorrect number of arguments entered!
Usage: ./foldXRun.py seqsDir mutsDir nativePDB pathToFoldX outputDir\n"""
        sys.exit(0)



