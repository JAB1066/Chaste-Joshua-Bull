# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 12:02:12 2017

@author: bull
"""
import sys
import os
import shutil
import string
# Script to combine different sections of the same simulation together. Designed to take a "burnt in spheroid" simulation which has been saved, and combine that with a simulation built upon that archive

baseFolderLocation = '/mi/share/scratch/bull/ChasteStuff/chaste_test_output/BenchmarkTumour/NecrosisPersists_12HoursOnly/'#sys.argv[1]
firstSimulationResultsFolder = 'results_from_time_120/'#'CombinedResults/'#'results_from_time_0/'#sys.argv[2] # results_from_time_0
secondSimulationResultsFolder = 'results_from_time_240/'#sys.argv[3] # results_from_time_120
startTimeInHoursOfFirstSimulation = 120
startTimeInHoursOfSecondSimulation = 240
outputTimestepFrequency = 120
offset = (startTimeInHoursOfSecondSimulation-startTimeInHoursOfFirstSimulation)*outputTimestepFrequency

combinedResultsFolder = baseFolderLocation + "CombinedResults/"
if ~os.path.exists(combinedResultsFolder) == -1:
    os.mkdir(combinedResultsFolder)

for filename in os.listdir(baseFolderLocation+firstSimulationResultsFolder):
    destName = combinedResultsFolder + filename
    sourceName = baseFolderLocation + firstSimulationResultsFolder + filename
    shutil.copyfile(sourceName,destName)
    #print(filename)
    
for filename in os.listdir(baseFolderLocation+secondSimulationResultsFolder):
    if filename.split('.')[-1] == 'vtu':
        #print(oldFilename)
        oldFilename = filename.split('.')[0]
        newFilenameNumber = int(oldFilename.split('_')[-1]) + offset
        newFilename = oldFilename.split('_')[:-1]
        newFilename.append(str(newFilenameNumber))
        newFilename = "_".join(newFilename)
        newFilename = newFilename + '.vtu'
        #print(newFilename)
        destName = combinedResultsFolder + newFilename
        sourceName = baseFolderLocation + secondSimulationResultsFolder + filename
        shutil.copyfile(sourceName,destName)
        #print(newFilename)
    



