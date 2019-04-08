## SY 23/8/18
## Copies delta files to new directory if suffix is < 5.
## Use following command line to copy selected files to new directory:
## xarg -I % --arg-file=filename.txt cp ./source-directory/% ./destination-directory/
## xarg -I % --arg-file=deltalist.txt cp ./ws_deltas/% ./ws_deltas_density1/

import glob

alldeltas = glob.glob("*.fits.gz")
ndel = len(alldeltas)
i=1
fout = open('deltalist.txt', 'w')

for filename in alldeltas:
    print(i, ndel)
    i+=1

    fin = filename[:-9]
    k = int(filename[-9:-8])
    if k < 5:
        fout.writelines(filename + '\n')

fout.close()


