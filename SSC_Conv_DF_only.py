##Cluster SSC web conversion function (handles 1 sc at once). checks file name for spacecraft to process correctly
##this one does NOT save as CSV

def SSCWEBConvNoSave(cdf_file):
    
    if 'cluster1' in cdf_file:
        sc_name = 'C1'
        cluster1 = cdflib.CDF(cdf_file)
        datetimes_1 = cdflib.cdfepoch.encode(cluster1['Epoch'])
        datetimes_series_1 = pd.Series(datetimes_1) 
        datetimes_a_1 = pd.to_datetime(datetimes_series_1, errors = 'coerce')
        datetimes_a_1 = datetimes_a_1.to_frame(name="datetime")
        
        df = pd.DataFrame({'X_gse': cluster1['XYZ_GSE'][:,0], 'Y_gse': cluster1['XYZ_GSE'][:,1], 'Z_gse': cluster1['XYZ_GSE'][:,2]})
        df = df.join(datetimes_a_1)
        df = df.set_index('datetime')
        
    if 'cluster2' in cdf_file:
        sc_name = 'C2'
        cluster2 = cdflib.CDF(cdf_file)
        datetimes_2 = cdflib.cdfepoch.encode(cluster2['Epoch'])
        datetimes_series_2 = pd.Series(datetimes_2) 
        datetimes_a_2 = pd.to_datetime(datetimes_series_2, errors = 'coerce')
        datetimes_a_2 = datetimes_a_2.to_frame(name="datetime")
        
        df = pd.DataFrame({'X_gse': cluster2['XYZ_GSE'][:,0], 'Y_gse': cluster2['XYZ_GSE'][:,1], 'Z_gse': cluster2['XYZ_GSE'][:,2]})
        df = df.join(datetimes_a_2)
        df = df.set_index('datetime')
        
    if 'cluster3' in cdf_file:
        sc_name = 'C3'
        cluster3 = cdflib.CDF(cdf_file)
        datetimes_3 = cdflib.cdfepoch.encode(cluster3['Epoch'])
        datetimes_series_3 = pd.Series(datetimes_3) 
        datetimes_a_3 = pd.to_datetime(datetimes_series_3, errors = 'coerce')
        datetimes_a_3 = datetimes_a_3.to_frame(name="datetime")
        df = pd.DataFrame({'X_gse': cluster3['XYZ_GSE'][:,0], 'Y_gse': cluster3['XYZ_GSE'][:,1], 'Z_gse': cluster3['XYZ_GSE'][:,2]})
        df = df.join(datetimes_a_3)
        df = df.set_index('datetime')
    
    if 'cluster4' in cdf_file:
        sc_name = 'C4'
        cluster4 = cdflib.CDF(cdf_file)
        datetimes_4 = cdflib.cdfepoch.encode(cluster4['Epoch'])
        datetimes_series_4 = pd.Series(datetimes_4) 
        datetimes_a_4 = pd.to_datetime(datetimes_series_4, errors = 'coerce')
        datetimes_a_4 = datetimes_a_4.to_frame(name="datetime")
        
        df = pd.DataFrame({'X_gse': cluster4['XYZ_GSE'][:,0], 'Y_gse': cluster4['XYZ_GSE'][:,1], 'Z_gse': cluster4['XYZ_GSE'][:,2]})
        df = df.join(datetimes_a_4)
        df = df.set_index('datetime')

    return(df)
