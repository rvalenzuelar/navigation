#!/usr/bin/python

# Replace coordinates in cfradial files
# 
# See http://www.unidata.ucar.edu/software/netcdf/examples/programs/
#
# Raul Valenzuela
# April, 2015

import sys 

# define function
def replace_cfradial_coords( stdtape_filepath ):

	from netCDF4 import Dataset
	import glob
	import numpy as np 
	import pandas as pd	

	# open standard tape file for reading
	stdtape_file = Dataset(stdtape_filepath,'r') 

	# get stdtape timestamp
	base_time=stdtape_file.variables['base_time'][:]
	stdtape_secs=stdtape_file.variables['Time'][:]
	stdtape_timestamp=pd.to_datetime(stdtape_secs+base_time,unit='s')
	stdtape_lats=stdtape_file.variables['LAT'][:]
	stdtape_lons=stdtape_file.variables['LON'][:]

	# close the file
	stdtape_file.close()

	# pandas dataframe for standar tape
	d={'lats':stdtape_lats,'lons':stdtape_lons}
	df_stdtape=pd.DataFrame(data=d,index=stdtape_timestamp)

	# get a list cfradial files
	nclist = glob.glob('cfrad.*')
	nlist=len(nclist)
	print "Folder contains ",str(nlist)," cfradial files"

	for f in np.arange(nlist):
		# open cfradial file for reading and writing
		print 'Processing: '+nclist[f]
		cfrad_file = Dataset(nclist[f],'r+') 

		# get cfradial timestamp
		start_datetime_nparray=cfrad_file.variables['time_coverage_start'][0:20]
		strpattern=''
		start_datetime_str=strpattern.join(str(v) for v in start_datetime_nparray)
		time_format="%Y-%m-%dT%H:%M:%SZ"
		cfrad_start_datetime=pd.to_datetime(start_datetime_str,format=time_format)
		cfrad_time = cfrad_file.variables['time'][:]
		cfrad_secs=pd.to_timedelta(cfrad_time.astype(int),unit='s')
		cfrad_timestamp=cfrad_start_datetime+cfrad_secs

		# remove duplicated timestamps (str)
		unique_timestamp=cfrad_timestamp.drop_duplicates()
		nstamps=unique_timestamp.nunique()

		# cfradial coordinates
		cfrad_lats = cfrad_file.variables['latitude'][:]
		cfrad_lons = cfrad_file.variables['longitude'][:]

		# pandas dataframe for cfradial file
		d={'lats':cfrad_lats,'lons':cfrad_lons}
		df_cfrad=pd.DataFrame(data=d,index=cfrad_timestamp)
		df_cfrad_new=df_cfrad.copy()

		for t in np.arange(nstamps):
			timestamp=str(unique_timestamp[t])
			new_lats=df_stdtape[timestamp]['lats']
			new_lons=df_stdtape[timestamp]['lons']
			df_cfrad_new.loc[timestamp,'lats']=new_lats
			df_cfrad_new.loc[timestamp,'lons']=new_lons

		cfrad_file.variables['latitude'][:]=df_cfrad_new['lats'].values
		cfrad_file.variables['longitude'][:]=df_cfrad_new['lons'].values

		# print ''
		# print df_cfrad['lats'].values
		# print ''
		# print type(df_cfrad_new['lats'].values)
		# print ''
		# print df_cfrad_new['2001-01-25 18:44:51']['lats']
		# print ''

		# close the file.
		cfrad_file.close()

# call function
stdtape=sys.argv[1]
replace_cfradial_coords(stdtape)