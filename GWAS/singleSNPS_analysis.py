import numpy as np
import subprocess
import json
import os
from os.path import isfile, join
from collections import defaultdict
import matplotlib.pyplot as plt
import scipy 
from scipy.spatial.distance import euclidean 

def HUQC_SSNP(settings):
    # Computes quality control in the merged files 
    # 1. Creates temporary files with the rsIDs correpsonding to MHC CHR 6 
    # 2. Filters out the MHC CHR 5 and maf at 0.05 and geno and missingness at 0.1
    # Output: filtHU in ./Data

    # Perform quality control
    print("Performing quality control")
    subprocess.call(['bash', 'qualityControl_singleSNPs.sh'])
    print("Quality control successfully performed")

def HUPCA_SSNP(settings):
    # Computes PCA on HU master files with common SNPs betwen cases and controsl 
    # no MHC and maf at 0.05 and geno and missingness at 0.1

    if bool(settings['performPCA']) is True:
        # Perform PCA
        print("Performing PCA")
        subprocess.call(['bash', 'pca_singleSNPs.sh'])
        print("PCA successfully performed")
    else:
        print("PCA deactivated")

def plotPCA(settings, f):

    # Initialize
    pc1 = list()
    pc2 = list()
    badge = list()

    # Read file
    with open(f, 'r') as inFile:
        for row in inFile:
            row = row.split(' ')
            pc1.append(float(row[2]))
            pc2.append(float(row[3]))
            badge.append(int(row[1].split('_')[1]))

    # Plot
    pc1 = np.array(pc1) ; pc2 = np.array(pc2) ; badge =  np.array(badge)
    fig, ax = plt.subplots()
    col = list()
    count = 0 
    for b in np.unique(badge):
        idx = np.where(badge == b)
        col.append(np.random.rand(3,))
        ax.scatter(x = pc1[idx], y = pc2[idx], c = col[count], \
            label = b, alpha = 0.3, linewidths=1)   
        count += 1
    plt.title("Principal components")
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.legend()
    plt.show(block=False)

    # Plot PCs individually 

    fig1, ax1 = plt.subplots(nrows=3, ncols=7)
    fig1.subplots_adjust(hspace=0.5)
    count = 0
    for ax, b in zip(ax1.flatten(), np.unique(badge)):
        idx = np.where(badge == b)
        ax.scatter(x = pc1[idx], y = pc2[idx], c = col[count], \
            label = b, alpha = 1, linewidths = 1)
        ax.set(title = "Plate " + str(b), xlabel = 'PC1', ylabel = 'PC2')
        count += 1
    plt.show()

def addPhenotype(settings, f, sep):
    # Opens filtHU and adds the phenotype based on the list of 
    # controls and cases. If not in cases or controls, ID is added 
    # to list and written in exclude.txt
    # Output modHU.fam 

    # Read GWAS IDS
    controlsID = list() 
    with open(settings['file']['GWASIDsControls'], 'r') as inFile:
        for row in inFile:
            controlsID.append(row.split('"')[1])
    casesID = list()
    with open(settings['file']['GWASIDsCases'], 'r') as inFile:
        for row in inFile:
            casesID.append(row.split('"')[1])

    # Feed phenotype to fam file and create exclude list
    exclude = list()
    with open(settings['file']['modHU_SSNP.fam'], 'w') as outFile:
        with open(f + '.fam', 'r') as inFile:
            for row in inFile:
                row = row.split(sep)
                subjectID = row[0].split('_')[3].split('.')[0]
                if subjectID in controlsID:
                    row[-1] = '1'
                    outFile.write('\t'.join(row) + '\n')
                elif subjectID in casesID:
                    row[-1] = '2'
                    outFile.write('\t'.join(row) + '\n')
                else:
                    exclude.append(row[0])
                    outFile.write('\t'.join(row))


def patientMatching(settings):

    # Perform patient matching
    print("Computing patient matching as the Euclidean distance between each case and control")

    # GET IDS and PCs
    patIDs = np.array([])
    patIDsGWAS = np.array([])
    PCs = list()
    IDs = list()
    with open(settings['file']['PCAeigenvectors_SSNP'], 'r') as inFile:
        for row in inFile:
            patIDsGWAS = np.append(patIDsGWAS, row.split(' ')[0])
            ID = np.asanyarray(row.split(' ')[0].split('_')[-1].split('.')[0])
            IDs.append(row.split(' ')[0].split('_')[-1].split('.')[0])
            patIDs = np.append(patIDs, ID)
            PCs.append(np.asanyarray([float(pc) for pc in row.split(' ')[2:]]))

    # Open cases ID
    casesID = np.array([])
    with open(settings['file']['GWASIDsCases'], 'r') as inFile:
        for row in inFile:
            casesID = np.append(casesID, row.split("\"")[1])

    # Open controls ID and get the ones in PatIDs (some removed in filtering)
    controlsID = np.array([])
    with open(settings['file']['GWASIDsControls'], 'r') as inFile:
        for row in inFile:
            controlsID = np.append(controlsID, row.split("\"")[1])
    controlsID = np.asanyarray([id for id in controlsID if id in patIDs])

    # Get matching patients
    ratio = settings['ControlCaseRatio']

    # For each case, compute euclidean distance
    keepIndxs = list()
    for case in casesID:

        # If case filtered or in exclusion list, continue 
        if case not in patIDs:
            continue
        
        # Get the index and PCs of the case. Store.
        idx = int(np.where(patIDs == case)[0][0])
        keepIndxs.append(idx)
        patPCs = PCs[idx]

        # Compute euclidean distance against all controls 
        patEuDist = list()
        for control in controlsID: 

            # If control filtered or in exclusion list, continue
            if control not in patIDs:
                continue

            # Get controls PCs and compute euclidean distance
            ctrlidx = int(np.where(patIDs == control)[0][0])
            ctrlPCs = PCs[ctrlidx]
            patEuDist.append(euclidean(patPCs, ctrlPCs))

        # Get lower distances between case and controls
        tempidx = np.asanyarray(patEuDist).argsort()[:ratio]
        ctrlsIDtemp = controlsID[tempidx]

        # Keep indxs
        for ctrl in ctrlsIDtemp:
            keepIndxs.append(np.where(patIDs == ctrl)[0][0])
        
    keepIndxs = np.unique(keepIndxs)
    # Write list of patients 
    patIDsGWAS = patIDsGWAS[keepIndxs]
    with open(settings['file']['patList_SSNP'], 'w') as outFile:
        for pat in patIDsGWAS:
            outFile.write(pat + ' ' + pat + '\n')

def postMatchPCA(settings):

    # Compute PCA filtering out non-matched patients
    print("Performing PCA after matching")
    subprocess.call(['bash','postPCA_singleSNPs.sh'])
    print("PCA finished")

    # Read PCs
    PCs = list(); IDs = list()
    with open(settings['file']['PCAeigenvectors_SSNP'],'r') as inFile:
        for row in inFile:
            PCID = [row.split(' ')[0], row.split(' ')[0]]
            PCID = PCID + [float(pc) for pc in row.split(' ')[2:]]
            PCs.append(PCID)
    with open(settings['file']['patList_SSNP'], 'r') as inFile:
        for row in inFile:
            IDs.append(row.split(' ')[0])
    PCID = [pc for pc in PCs if pc[0] in IDs]

    with open(settings['file']['postPCAeigenvectors_SSNP'], 'w') as outFile:
        for pc in PCID:
            pcSTR = ' '.join(pc[0:2])
            for val in pc[2:-1]: pcSTR = pcSTR + ' {:}'.format(val)
            outFile.write(pcSTR + '\n')
        
    f = settings['file']['postPCAeigenvectors_SSNP']
    plotpostPCA(settings, f)

def plotpostPCA(settings, f):

    # Initialize
    pc1 = list()
    pc2 = list()
    IDs = list()

    # Read file
    with open(f, 'r') as inFile:
        for row in inFile:
            row = row.split(' ')
            pc1.append(float(row[2]))
            pc2.append(float(row[3]))
            IDs.append(str(row[1].split('_')[-1].split('.')[0]))
    
    # Load cases and controls ID
    controlsID = list() 
    with open(settings['file']['GWASIDsControls'], 'r') as inFile:
        for row in inFile:
            controlsID.append(row.split('"')[1])
    casesID = list()
    with open(settings['file']['GWASIDsCases'], 'r') as inFile:
        for row in inFile:
            casesID.append(row.split('"')[1])

    # Get label (case vs control)
    label = list()
    for ID in IDs:
        if ID == '86A10':
             print('stop')
        if ID in casesID:
            label.append('cases')
        elif ID in controlsID:
            label.append('controls')
        else:
            label.append('exclude')
    if 'exclude' in label:
        label = np.array(label) ; excludeIdx = np.where(label == 'exclude'); label = np.delete(label, excludeIdx)
        pc1 = np.delete(np.array(pc1), excludeIdx) ; pc2 = np.delete(np.array(pc2), excludeIdx)
    else:
        pc1 = np.array(pc1) ; pc2 = np.array(pc2) ; label =  np.array(label)

    # Plot
    fig, ax = plt.subplots()
    col = list()
    count = 0 
    for b, c in zip(np.unique(label), ["#33CC33", "#0066CC"]):
        idx = np.where(label == b)
        col.append(np.random.rand(3,))
        ax.scatter(x = pc1[idx], y = pc2[idx], c = c, \
            label = b, alpha = 0.7, linewidths=1)   
        count += 1
    plt.title("Principal components")
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.legend()
    plt.show(block=False)

    # Plot PCs individually 

    fig1, ax1 = plt.subplots(nrows=1, ncols=2)
    fig1.subplots_adjust(hspace=0.5)
    count = 0
    for ax, b, c in zip(ax1.flatten(), np.unique(label), ["#33CC33", "#0066CC"]):
        idx = np.where(label == b)
        ax.scatter(x = pc1[idx], y = pc2[idx], c = c, \
            label = b, alpha = 1, linewidths = 1)
        ax.set(title = "Plate " + str(b), xlabel = 'PC1', ylabel = 'PC2')
        count += 1
    plt.show()

def filtPostMatchData(settings):

    # Initialize
    PCs = list()
    IDs = list()
    f = settings['file']['postPCAeigenvectors_SSNP']

    # Import PCs
    with open(f, 'r') as inFile:
        for row in inFile:
            row = row.split(' ')
            PCs.append([float(pc) for pc in row[2:]])
            IDs.append(row[1:2])

    # Filter PCs using thresholds
    idxs = list()
    count = 0

    if settings['diagnosis'] == 'LGI1':
        for PC in PCs:
            if PC[1] > 0.0229 and PC[1] < 0.057:
                count += 1
                continue
            elif PC[0] > 0.0014 and PC[1] < 0.04:
                count += 1
                continue
            else:
                idxs.append(count)
                count += 1
    elif settings['diagnosis'] == 'HU':
         for PC in PCs:
            if PC[1] < 0.02534 and PC[0] > -0.00895:
                count +=1 
                continue
            if PC[1] > 0.04 and PC[1] < 0.08: 
                count += 1 
                continue
            else:
                idxs.append(count)
                count += 1
    elif settings['diagnosis'] == 'Caspr':
        for PC in PCs:
            if PC[0] > -0.0037:
                count += 1
                continue
            else:
                idxs.append(count)
                count += 1
    elif settings['diagnosis'] == 'YO':
        for PC in PCs:
            if PC[0] > 0.1:
                count += 1
                continue
            else:
                idxs.append(count)
                count += 1
    
    # Rewrite PCA 
    PCs = list(list(PC) for PC in np.asanyarray(PCs)[idxs]); 
    IDs = list(ID[0] for ID in np.asanyarray(IDs)[idxs])
    with open(settings['file']['filtposPCAeigenvectors_SSNP'], 'w') as outFile:
        for PC, ID in zip(PCs, IDs):
            PC = [str(x) for x in PC]
            outFile.write(ID + ' ' + ID + ' ' + ' '.join(PC) + '\n')

    # Get list of patients and write 
    with open(settings['file']['filtPostMatchPatList_SSNP'], 'w') as outFile:
        for pat in IDs:
            outFile.write(pat + ' ' + pat + '\n')

    # Plot PCA 
    f = settings['file']['filtposPCAeigenvectors_SSNP']
    plotpostPCA(settings, f)


def associatioAnalysis(settings):

    # Perform association analysis
    print("Performing Asssociation analysis")
    subprocess.call(['bash','associationAnalysis_SSNP.sh'])
    print ('\n')
    print("Association analysis performed successfully")


def main():

    # Open settings
    with open('settings.json', 'r') as inFile:
        settings = json.load(inFile)

    # Compute quality control
    HUQC_SSNP(settings)

    # Add phenotype to filtered 
    f = settings['file']['filtHU_SSNP']
    sep = ' '
    addPhenotype(settings, f, sep)

    # Compute PCA 
    HUPCA_SSNP(settings)

    # Plot PCA
    f = settings['file']['PCAeigenvectors_SSNP']
    plotpostPCA(settings, f)

    # Patient matching
    patientMatching(settings)

    # Plot PCA
    postMatchPCA(settings)

    # Manually remove postMatch data
    filtPostMatchData(settings)

    # Add field separator in PCA
    # TODO: CHANGE THIS
    with open(settings['file']['PCAeigenvectors_SSNP'], 'r') as inFile: 
        with open('./ResourcesGWAS/PCAtst_SSNP.eigenvec', 'w') as outFile:
            header = ['PC' + str(x) for x in range(1,21)]
            header = ['FID','IID'] + header
            outFile.write('\t'.join(header) + '\n') 
            for row in inFile:
                row = row.split(' ')
                row = '\t'.join(row)
                outFile.write(row)

    # Perform association analysis
    associatioAnalysis(settings)

    # Manhattan Plot 
    subprocess.run(['Rscript', '--vanilla','manhattanPlot.R','singleSNP'])

    
if __name__ == "__main__":
    main()    
    
