#extract tarfiles
import tarfile
import glob

def tarf_extract(tarfilelist):
    error_list = []

    for i in tarfilelist:
        try:
            with tarfile.open(i) as tar:
                tarname = tar.getnames()
                tar.extractall(path='/data/scratch/apx059/23_Years_Data/CDFs')
        except:
            error_list.append(i)
            continue

##unzip Cluster tgzs

list_all = []

path = "/data/scratch/apx059/23_Years_Data/C1/**"

for path in glob.glob(path, recursive=True):
    list_all.append(path)
    
#list with only files, not folders
tgz_file_list = []

for element in list_all:
    if '.tgz' in element:
        tgz_file_list.append(element)
        
tarf_extract(tgz_file_list)

##unzip Cluster tgzs

list_all = []

path = "/data/scratch/apx059/23_Years_Data/C2/**"

for path in glob.glob(path, recursive=True):
    list_all.append(path)
    
#list with only files, not folders
tgz_file_list = []

for element in list_all:
    if '.tgz' in element:
        tgz_file_list.append(element)
        
tarf_extract(tgz_file_list)

##unzip Cluster tgzs

list_all = []

path = "/data/scratch/apx059/23_Years_Data/C4/**"

for path in glob.glob(path, recursive=True):
    list_all.append(path)
    
#list with only files, not folders
tgz_file_list = []

for element in list_all:
    if '.tgz' in element:
        tgz_file_list.append(element)
        
e_list = tarf_extract(tgz_file_list)
print(e_list)