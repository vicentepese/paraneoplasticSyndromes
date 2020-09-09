import numpy as np
import json
import sys
import subprocess

def getrsIDs(settings):

    # Read and get IDs controls
    controls = list()

    with open(settings['file']['HU_controls.map'], 'r') as inFile:
        print('Getting IDs from controls')
        for i in range(0,4): next(inFile)   # Skip headers
        for row in inFile:
            row = row.split('\t')
            controls.append(row[1])
        print('Done')

    # Read and get IDs cases
    cases = list()
    with open(settings['file']['HU_cases.map'], 'r') as inFile:
        print('Getting IDs from cases')
        for i in range(0,4): next(inFile)   # Skip headers
        for row in inFile:
            row = row.split('\t')
            cases.append(row[1])
        print('Done')
    
    # Find common ids
    cmnrsIDs = list(set(cases) & set(controls))

    return controls, cases

def filterSNPs(settings, controls, cases):

    # Find indexes of commmon rsIDs in cases: 
    notCmnrsIDs = list(set(controls) -  set(cases))

    # Write 
    with open(settings['file']['notComnrsIDs.txt'], 'w') as outFile:
        for rsid in notCmnrsIDs:
            outFile.write(rsid + '\n')

    # Run plink 
    print('Filtering SNPS')
    result = subprocess.run(['bash', './Utils/plink_filterSNPs.sh'], stdout=subprocess.PIPE)
    print('Done')
    
def main():

    # Load options
    with open("settings.json", 'r') as inFile:
        settings = json.load(inFile)

    # Get common rsIDs from patients and controls
    controls, cases = getrsIDs(settings)
    
    # Print 
    print("Length of common rsIDs is equal to the length of cases")
    subprocess.call(['cp', settings['file']['HU_cases.map'], '../../HU.map'], stdout=subprocess.PIPE)
    print("HU.map is HU_cases.map")

    # Filter SNPS in cases:
    filterSNPs(settings, controls, cases)


    pass

if __name__ == "__main__":
    main()