#extract tarfiles
import tarfile
import glob

def tarf_extract(tarfilelist):
    error_list = []

    for i in tarfilelist:
        try:
            with tarfile.open(i) as tar:
                tarname = tar.getnames()
                tar.extractall(path='/data/SPCS-HIETALA-Shocks/GIPM-MAPPING/Cluster_CDFs/')
        except:
            error_list.append(i)
            continue
            
    print(error_list)  
          

list_all = []

#folder with input tarfiles

path = "/data/SPCS-HIETALA-Shocks/GIPM-MAPPING/Cluster_Source_Tarfiles/**"

for path in glob.glob(path, recursive=True):
    list_all.append(path)
    
#ensure only trying to unpack tgz files

tgz_file_list = []

for element in list_all:
    if '.tgz' in element:
        tgz_file_list.append(element)
        
tarf_extract(tgz_file_list)

