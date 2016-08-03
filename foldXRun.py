#!/usr/bin/env python
from fabric.api import local, env, settings, lcd, hide
import re
from tempfile import mkdtemp
from shutil import rmtree, copy
from os.path import isfile, join, normpath, exists, abspath, dirname, basename, splitext
import sys
from os import listdir, makedirs
import numpy
import csv

repairFile = """\
command=RepairPDB
pdb=%s
"""
stabilizeFile = """\
command=Stability
pdb=%s
"""

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


    # Setup repair PDBs
    repairJobPDBs = []
    for modelFile in seqModels:
        repairJobPDBs.append("seq_%s" % modelFile)
    for modelFile in mutModels:
        repairJobPDBs.append("mut_%s" % modelFile)
    repairJobPDBs.append("native.pdb")
    with open(join(scratchDir, "repairPDBs"), 'w') as f:
        f.write("\n".join(repairJobPDBs))


    # Now we handle the stability job in the same way.


    # Setup stabilize PDBs
    stabilizeJobPDBs = []
    for modelFile in seqModels:
        stabilizeJobPDBs.append("seq_%s_Repair.pdb" % splitext(modelFile)[0])
    for modelFile in mutModels:
        stabilizeJobPDBs.append("mut_%s_Repair.pdb" % splitext(modelFile)[0])
    stabilizeJobPDBs.append("native_Repair.pdb")
    with open(join(scratchDir, "stabilizePDBs"), 'w') as f:
        f.write("\n".join(stabilizeJobPDBs))

    # Setup report file list
    reportSuffix = "_0_ST.fxout"
    reportFiles = map(lambda pdb: splitext(pdb)[0] + reportSuffix,  stabilizeJobPDBs)
    


    repairRunFileText = repairFile
    stabilizeRunFileText = stabilizeFile
    with open(join(scratchDir, 'repairRunFile.cfg'), 'w') as repair:
        with open(join(scratchDir, 'stabilizeRunFile.cfg'), 'w') as stabilize:
            repair.write(repairRunFileText)
            stabilize.write(stabilizeRunFileText)
    with lcd(scratchDir):
        with hide('everything'):
            print "Running repair job."
            print repairJobPDBs
            for repairJobPDB in repairJobPDBs:
                local("%s --command=RepairPDB --pdb=%s" % (foldXPath, repairJobPDB))
            print "Running stabilize job."
            for stabilizeJobPDB in stabilizeJobPDBs:
                local("%s --command=Stability --pdb=%s" % (foldXPath, stabilizeJobPDB))
            print "Done~"

    print results
    print "FoldX Complete - Now running analysis on output"
    reportPaths = map(lambda reportFile: join(dataDir, reportFile), reportFiles)
    reportAnalyze(scratchDir, join(dataDir, 'report.txt'), basename(normpath(seqDirectory)))

def reportAnalyze(reportPath, outputPath, reportName):
    # From  http://foldxsuite.crg.eu/command/Stability#inout-parameters
    fieldHeaders = ["name", "total energy", "Backbone Hbond", "Sidechain Hbond", "Van der Waals", "Electrostatics", "Solvation Polar", "Solvation Hydrophobic", "Van der Waals clashes", \
               "Entropy Side Chain", "Entropy Main Chain", "Sloop Entropy", "Mloop Entropy", "Cis Bond", "Torsional Clash", "Backbone Clash", "Helix Dipole", "Water Bridge", \
               "Disulfide", "Electrostatic Kon", "Partial Covalent Bonds", "Energy Ionisation", "Entropy Complex", "Residue Number"]

    native = None
    seqs = []
    muts = []
    seqsSTDs = {}
    mutsSTDs = {}
    seqsMeans = {}
    mutsMeans = {}

    
    reportFiles = [f for f in listdir(reportPath) if re.search(".*_ST.fxout", f) and isfile(join(reportPath, f))]

    # Reading in values from the Stabilized report
    for reportFile in reportFiles:
        with open(join(reportPath, reportFile), 'r') as f:
            reader = csv.DictReader(f, delimiter='\t', fieldnames=fieldHeaders)
            row = reader.next()
            if re.match(".*mut", row["name"]):
               muts.append(row)
            elif re.match(".*seq", row["name"]):
                seqs.append(row)
            else:
                native = row

    # Python is really bad!!
#    for i in seqs:
#        i['name'] = i['']
#        del i['']
#    for i in muts:
#        i['name'] = i['']
#        del i['']
#    native['name'] = native['']
#    del native['']

    # Calculate means, stds for each value
    for key in seqs[0]:
        # We don't want stats on these!
        if key in ('name', 'Residue Number'):
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



