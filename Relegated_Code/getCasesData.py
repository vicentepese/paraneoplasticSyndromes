import numpy as np
import os 
from collections import defaultdict
import csv 
import json 
import subprocess


def getIds(settings):

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
    Dx = paraNeoSample['Dx']
    ids = list()
    idsOut = list()

    for ID, dx in zip(FIDs, Dx):
        if dx == 'HU':
            if len(ID) == 5:
                ids.append(ID[0:2])
            else:
                ids.append(ID[0:3])
            idsOut.append(ID)

    # Get unique IDs
    uniqueIDs = [int(ID) for ID in np.unique(ids)]
    uniqueIDs.sort()
    uniqueIDs = ''.join([str(ID) + ', ' for ID in uniqueIDs])

    # Print 
    print('GWAS are located in plates:' + uniqueIDs)

    # Return 
    return idsOut

def prepCases(settings, ids):

    # Get file names 
    files = os.listdir(settings['directory']['GWAS_cases']) 

    if settings['file']['HU_cases.ped'].split('/')[-1] not in files:
        print('Preparing cases data')
        count = 0
        for file in files:
            if '.ped' in file:
                print('Current file: ' + file)
                with open(settings['file']['HU_cases.ped'], 'a') as outFile:
                    with open(settings['directory']['GWAS_cases'] + file, 'r') as inFile:
                        for row in inFile:
                            row = row.split('\t')
                            patID = row[0].split('_')[-1].split('.')[0]
                            if len(patID) > 2 and patID in ids:
                                if 'Copy' not in row[1]:
                                    row[0] = patID
                                    row.insert(1, str(2))
                                    row = '\t'.join(row)
                                    outFile.write(row)
                                    count += 1

        
        # Verbose
        print('Cases data preparation finished')
        print('\n')

    else: 
        count = subprocess.run(['bash', settings['util']['countRows'],\
             settings['file']['HU_cases.ped']], stdout=subprocess.PIPE)
        count = int(count.stdout.decode('utf-8').split('\n')[0])
        
    print('Number of patients: ' + str(count))

def main():

    # Import settings 
    with open('settings.json', 'r') as inFile:
        settings = json.load(inFile)

    # Get plates
    ids = getIds(settings)

    # Get cases data
    prepCases(settings, ids)
    


if __name__ == "__main__":
    main()