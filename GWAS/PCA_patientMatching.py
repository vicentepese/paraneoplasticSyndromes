import numpy as np
import json 
from collections import defaultdict, Counter
import csv
import matplotlib.pyplot as plt
from scipy.spatial import distance

########## DESCRIPTION #########
# Case-control Principal Components matching through euclidean distance 
# The data provided is: 
#   - "GWAS/ResourcesGWAS/HLA_CALLS_LGL1_CASE_CONTROLS.csv": A datasheet with the HLA calls of 1353 subjects. Contains diagnosis (case-control) and IDs. 
#       The file is used as a reference. 
#   - "GWAS/ResourcesGWAS/GPC_PMRA_upto121.eigenvec": A dataset with the PCs of 9001 subjects, including the 1353 subjects aforementioned.
#   - "GWAS/ResourcesGWAS/GPC_PMRA_upto121_PIHAT_025.genome": IBD of the 9001 aforementioned subjects. 
################################

def preprocessData(settings):

    # Initialize HLA calls, PCs and genome dicts. IDs in GWAS format
    HLA_calls, PCs, genomeData = defaultdict(lambda: defaultdict(str)), \
         defaultdict(lambda: defaultdict(str)), list()

    # Import HLA calls 
    with open(settings['file']['HLA_calls_all_plates'],'r') as inFile:
        csv_reader = csv.DictReader(inFile, delimiter = ',')
        for row in csv_reader:
            HLA_calls[row['GWASID']] = row; HLA_calls[row['GWASID']].pop('GWASID', None)

    # Extract cases 
    HLA_calls_cases = {ID: HLA_calls[ID] for ID in HLA_calls if int(HLA_calls[ID]['Dx']) == 1}

    # Import PCs
    with open(settings['file']['PCs_all_plates'], 'r') as inFile:
        csv_reader = csv.DictReader(inFile,  delimiter = ' ')
        next(csv_reader)
        for row in csv_reader:
            if 'Stanford' in row['IID']:
                IID = row['IID'].split('.')[0].split('_')[-1]
            else:
                IID = row['IID']
            PCs[IID] = row
            PCs[IID].pop('IID', None)

    # Import genome 
    with open(settings['file']['genome_all_plates'],'r') as inFile:
        csv_reader = csv.DictReader(inFile, delimiter = ',')
        for row in csv_reader:
            if 'Stanford' in row['IID1']:
                IID = row['IID1'].split('.')[0].split('_')[-1]
            else:
                IID = row['IID1']
            row.pop('IID1', None)
            genomeData.append({IID: row }) 

    # Compute intersection by subject between HLA calls and PCs, and the IBD 
    commonIDs = list(set(list(HLA_calls.keys())).intersection(set(list(PCs.keys()))))
    HLA_calls_subset = {ID: HLA_calls[ID] for ID in HLA_calls.keys() if ID in commonIDs}
    PCs_subset = {ID: PCs[ID] for ID in PCs.keys() if ID in commonIDs}
    genomeData_subset = [subj for subj in genomeData if list(subj.keys())[0] in commonIDs]
    
    # Check IBD and remove duplicates 
    PI_HAT = [float(subj[list(subj.keys())[0]]['PI_HAT']) for subj in genomeData_subset]
    subj_dup = list()
    for subject in genomeData_subset:
        if float(subject[list(subject.keys())[0]]['PI_HAT']) >= settings['IBD_threshold']:
            subj_dup.append({"IID1": list(subject.keys())[0], \
                             "IID2": subject[list(subject.keys())[0]]['IID2'], \
                             "PI_HAT": subject[list(subject.keys())[0]]['PI_HAT']})
    
    # Count number of duplicates 
    subj_count = defaultdict(int)
    for subject in subj_dup:
        IID1, IID2 = subject['IID1'], subject['IID2']
        if 'Stanford' in IID1:
            IID1 = IID1.split('.')[0].split('_')[-1]
        if 'Stanford' in IID2: 
            IID2 = IID2.split('.')[0].split('_')[-1]
        subj_count[IID1] += 1
        subj_count[IID2] += 1
    count = np.array([subj_count[subj] for subj in subj_count])

    subj2remove_list = list()
    while len(np.unique(count)) !=1 and max(count) != 1:

        # Get max subject and append to list to remove 
        subj2remove = max(subj_count, key=subj_count.get)
        subj2remove_list.append(subj2remove)

        # Re-count number of duplicates
        remIdx = list()
        for (idx,pairSubjs) in enumerate(subj_dup):
            if pairSubjs['IID1'] == subj2remove or pairSubjs['IID2'] == subj2remove:
                remIdx.append(idx)
        for idx in sorted(remIdx, reverse=True): del subj_dup[idx]
            
        # Re-count 
        subj_count = defaultdict(int)
        for subject in subj_dup:
            IID1, IID2 = subject['IID1'], subject['IID2']
            if 'Stanford' in IID1:
                IID1 = IID1.split('.')[0].split('_')[-1]
            if 'Stanford' in IID2: 
                IID2 = IID2.split('.')[0].split('_')[-1]
            subj_count[IID1] += 1
            subj_count[IID2] += 1
        count = np.array([subj_count[subj] for subj in subj_count])

    # Remove half of the remaining subjects and add to list
    subj2remove_list += [ID['IID1'] for ID in subj_dup] 

    # Remove data   
    for subj in subj2remove_list:
        HLA_calls_subset.pop(subj, 'None')
        PCs_subset.pop(subj, 'None')

    # Return 
    return(HLA_calls_subset, PCs_subset, genomeData_subset)

def PCMatching(settings, HLA_calls, PCs, genomeData):

    # GWAS IDs of cases and controls 
    GWASID_cases, GWASID_controls = list(), list()
    for subj in HLA_calls:
        if int(HLA_calls[subj]['Dx'])== 1:
            GWASID_cases.append(subj)
        elif int(HLA_calls[subj]['Dx']) == 0:
            GWASID_controls.append(subj)
    GWASID_controls = np.asarray(GWASID_controls)

    # Initialize 
    controlsMatched = []

    # For each case
    for case in GWASID_cases:

        # Get PCs (1 to 2)
        casePCs = np.array([float(PCs[case]['PC1']), float(PCs[case]['PC2'])])

        # Compute euclidean distance with other controls:
        caseControlDists = list()
        for control in GWASID_controls:
            controlPCs = np.array([float(PCs[control]['PC1']), float(PCs[control]['PC2'])])
            eucl_dist = distance.euclidean(casePCs, controlPCs)
            caseControlDists.append(eucl_dist)

        # Get the closest number of controls as defined by settings 
        controlsMatched += GWASID_controls[np.argsort(caseControlDists)[0:settings['ControlCaseRatio']]].tolist()

    # Get unique controls 
    controlsMatched = np.unique(np.asanyarray(controlsMatched))

    # Return 
    return (controlsMatched)

def plotPCA(settings, controlsMatched, PCs, HLA_calls):

    # Subset data 
    PC1, PC2, label = list(), list(), list()
    patList = controlsMatched.tolist() + [ID for ID in HLA_calls if int(HLA_calls[ID]['Dx']) == 1]
    for pat in patList:
        PC1, PC2 = PC1 + [float(PCs[pat]['PC1'])], PC2 + [float(PCs[pat]['PC2'])]
        label += ['control' if pat in controlsMatched else 'case']

    # Count cases 
    num_cases =Counter(label)['case']
    num_controls = len(label) - num_cases

    # Plot
    fig, ax = plt.subplots()
    count = 0 
    for b, c in zip(['control', 'case'], ["#0066CC", "#33CC33"]):
        idx = np.where(np.asarray(label) == b)
        ax.scatter(x = np.asarray(PC1)[idx], y = np.asarray(PC2)[idx], c = c, \
            label = b, alpha = 0.7, linewidths=1)   
        count += 1
    plt.title("Principal components // {:} cases vs. {:} controls".format(num_cases, num_controls))
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.legend()
    plt.show()

     # Write controls Matched
    with open(settings['file']['patList_allplates'], 'w') as outFile:
        for pat in patList:
            outFile.write(pat + '\n')

def main():

    # Import settings 
    with open('GWAS/settingsGWAS.json','r') as inFile:
        settings = json.load(inFile)

    # Pre-process data 
    HLA_calls, PCs, genomeData = preprocessData(settings)

    # Compute case-control matching 
    controlsMatched = PCMatching(settings, HLA_calls, PCs, genomeData)

    # Plot PCA-matched cases and controls 
    plotPCA(settings, controlsMatched, PCs, HLA_calls)

if __name__ == "__main__":
    main()
