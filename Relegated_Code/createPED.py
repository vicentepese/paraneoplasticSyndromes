import json
import csv
import numpy as np 
from collections import defaultdict
import os 

def getPEDFields(settings):

    # Initialize inFile dictionnary
    paraNeoSample = defaultdict(list)

    # Load data 
    with open(settings['file']['paraneoplastic'], 'r') as inFile:
        for row in inFile: 
            if paraNeoSample:
                for value, var in zip(row.split(';'), varlist):
                    paraNeoSample[var].append(value)
            else:
                varlist = row.split(';')
                for var in varlist: paraNeoSample[var] = []
   
    # Initialize variables and output 
    FIDs = paraNeoSample['GWAS ID']
    Sex = paraNeoSample['Sex']
    Dx = paraNeoSample['Dx']
    ids = list()
    phenoInfo = defaultdict(list)

    for ID, dx, sMF in zip(FIDs, Dx, Sex):
        if dx == 'HU':
            if len(ID) == 5:
                ids.append(ID[0:2])
            else:
                ids.append(ID[0:3])
            phenoInfo['FID'].append(ID)
            phenoInfo['IID'].append(0)
            phenoInfo['IIDF'].append(0)
            phenoInfo['IIDM'].append(0)
            if sMF == 'M':
                phenoInfo['sex'].append(1)
            elif sMF == 'F':
                phenoInfo['sex'].append(2)
            else:
                phenoInfo['sex'].append(0)
            phenoInfo['phenotype'].append(1)

    # Get unique IDs
    uniqueIDs = [int(ID) for ID in np.unique(ids)]
    uniqueIDs.sort()
    uniqueIDsOut = uniqueIDs
    uniqueIDs = ''.join([str(ID) + ', ' for ID in uniqueIDs])

    # Print 
    print('GWAS are located in plates:' + uniqueIDs)

    return phenoInfo['FID'], phenoInfo

def createPEDCases(settings, idsOut, phenoInfo):

    # Get file names 
    files = os.listdir(settings['directory']['GWAS_cases'])
    
    # If .ped files not pre-processed
    if not any('.ped' in ff for ff in [f for f in os.listdir(settings['directory']['GWAS'])]):
        print('Pre-processing .ped files')
        # For each file, check patient ID if ped file and write output
        count = 0
        for file in files:
            if '.ped' in file:
                print('Current file: ' + file)
                with open(settings['file']['HU_cases.ped'], 'a') as outFile:
                    with open(settings['directory']['GWAS_cases'] + file, 'r') as inFile:
                        for row in inFile:
                            patID = row.split('\t')[0]
                            patID = patID.split('_')
                            if len(patID) > 2:
                                patID = patID[1] + patID[2]
                                if patID in idsOut:
                                    row = row.split('\t')
                                    outrow = [phenoInfo[key][count] for key in phenoInfo.keys()]
                                    row[0] = outrow[0]
                                    for i in range(1, len(outrow)):
                                        row.insert(i, str(outrow[i]))
                                    row = '\t'.join(row)
                                    outFile.write(row)
                                    count += 1
        # Verbose
        print('Number of patients: ' + str(count))
        print('Pre-processing finished')
        print('\n')

   

def mergeMap(settings):

    # Get file names 
    files = os.listdir(settings['directory']['GWASOG'])

    # Check 
    if not settings['file']['HU.map'] in files:
        print("Merging map files")
        with open(settings['file']['HU.map'], 'a') as outFile:
            for file in files:
                if '.map' in file:
                    print("Current file: " + file)
                    with open(settings['directory']['GWAS_cases'] + file, 'r') as inFile:
                        for row in inFile:
                            outFile.write(row)
    
        # Verbose
        print('Merging finished')

def main():
    
    # Load settings
    with open('settings.json', 'r') as inFile:
        settings = json.load(inFile)

    # Get plates
    ids, phenoInfo = getPEDFields(settings)

    # Get IDs
    createPEDCases(settings, ids, phenoInfo)

    # Merge map files 
    mergeMap(settings)





if __name__ == "__main__":
    main()