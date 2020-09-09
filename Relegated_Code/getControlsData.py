import numpy as np
import os 
from collections import defaultdict
import csv 
import json 
import subprocess

def prepControls(settings):

    # Load GWASIDs
    with open(settings['file']['GWASIDs'], 'r') as inFile:
        GWASIDs = list()
        for row in inFile:
            GWASIDs.append(row.split("\"")[1])

    # Get files names 
    files = os.listdir(settings['directory']['GWAS_OG']) 

    # Get patients ID
    patIDs = list()

    # Create controls data
    if settings['file']['HU_controls.ped'].split('/')[-1] not in files:
        print('Preparing cases data')
        count = 0
        for file in files:
            if '.ped' in file:
                print ('Current file: ' + file)
                with open(settings['file']['HU_controls.ped'], 'a') as outFile:
                    with open(settings['directory']['GWAS_OG'] + file, 'r') as inFile:
                        for row in inFile:
                            row = row.split('\t')
                            patID = row[0].split('_')[-1].split('.')[0]
                            if patID in GWASIDs and 'Copy' not in row[1] and patID not in patIDs:
                                row[0] = patID
                                patIDs.append(patID)
                                row.insert(1, str(1))
                                row = '\t'.join(row)
                                print(row[0:20])
                                outFile.write(row)
                                count += 1

        # Verbose
        print('Controls data preparation finished')
        print('\n')

    else:
        count = subprocess.run(['bash', settings['util']['countRows'],\
             settings['file']['HU_controls.ped']], stdout=subprocess.PIPE)
        count = int(count.stdout.decode('utf-8').split('\n')[0])

    print('Number of controls: ' + str(count))


    for file in files:
        if '.map' in file:
            count = subprocess.run(['bash', settings['util']['countRows'],\
                settings['file']['HU_controls.ped']], stdout=subprocess.PIPE)
            count = int(count.stdout.decode('utf-8').split('\n')[0])
            print(file + ' -- Number of rows: ' + str(count))

def main():

    # Load settings 
    with open('settings.json', 'r') as inFile:
        settings = json.load(inFile)

    # Create PED file
    prepControls(settings)
    

if __name__ == "__main__":
    main()
