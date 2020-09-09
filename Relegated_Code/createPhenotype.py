import json 
import pandas as pd
import csv

def filterSample(settings):

    # Initialize output
    filtSamples = list()

    # Read and filter 
    with open(settings['file']['paraneoplastic'], 'r') as inFile:
        next(inFile)
        for row in inFile:
            vars = row.split(';')

            # Label diagnosis
            if vars[2] == settings['diagnosis']:
                dx = 1
            else:
                dx = 0

            # Append
            filtSamples.append([vars[0], 1, dx])
    
    # Write output
    with open(settings['file']['phenotypeForward'], 'w') as outFile:
        writer = csv.writer(outFile, delimiter='\t')
        writer.writerows(filtSamples)

def getPatientID(settings):

    # Initialize 
    tst = list()

    # Read
    with open(settings['file']['Forward'] + '.ped', 'r') as inFile:
        for row in inFile:
            vars = row.split('\t')
            tst.append(vars[1] + '\n')
    
    # Write 
    with open('tst.csv', 'w') as outFile:
        outFile.writelines(tst)

def main():

    # Read settings
    with open('settings.json','r') as inFile:
        settings = json.load(inFile)
    
    # Read and filter Samples
    filterSample(settings)

    # Get patientID
    getPatientID(settings)


if __name__ == "__main__":
    
    main()