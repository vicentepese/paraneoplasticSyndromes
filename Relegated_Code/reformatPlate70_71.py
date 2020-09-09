import json


def main():

    # Import settings
    with open('settings.json','r') as inFile:
        settings = json.load(inFile)

    # Get files
    path = settings['directory']['GWAS_casesOG']
    files = [path + 'Plate70_Narcolepsy_updated.ped' , \
        path + 'Plate71_Narcolepsy_updated.ped']
    files_out = [path + 'Plate70_Narcolepsy_updated_mod.ped', \
        path + 'Plate71_Narcolepsy_updated_mod.ped']

    # Change FID name format
    for fileIn, fileOut in zip(files, files_out): 
        with open(fileIn, 'r') as inFile:
            with open(fileOut, 'w') as outFile:
                for row in inFile:
                    row = row.split(' ')
                    del row[0]; del row[1:5]
                    FID = row[0].split('_')
                    FID = 'Stanford_' + FID[1] + '_' + FID[2] + '_' + FID[1] + FID[2] + '.CEL'
                    row[0] = FID
                    outFile.write(' '.join(row) + '\n')


if __name__ == "__main__":
    main()