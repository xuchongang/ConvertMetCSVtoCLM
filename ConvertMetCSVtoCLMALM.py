# =======================================================================================
# =======================================================================================

import numpy as np
import sys
import getopt
import code  # For development: code.interact(local=locals())
from datetime import datetime
from matplotlib.dates import date2num, num2date
import csv
from scipy.io import netcdf
import matplotlib.pyplot as plt

simple_verify_ts = True

class ctrltype: 
    
    def __init__(self,csv_file,time_res_sec,n_header_row,n_fields, \
                 nodata_flags,grid_name_out,fill_value,missing_value):

        self.csv_file      = csv_file
        self.time_res_sec  = time_res_sec
        self.n_header_row  = n_header_row
        self.n_fields      = n_fields
        self.nodata_flags  = nodata_flags
        self.grid_name_out = grid_name_out
        self.fill_value    = fill_value
        self.missing_value = missing_value


class contype:

    def __init__(self,name,long_name,units,mode,value):
        
        self.name      = name
        self.long_name = long_name
        self.units     = units
        self.mode      = mode
        self.value     = value

class vartype:

    def __init__(self,name,long_name,units,mode,col_id,unit_mult,unit_off):
    
        self.name      = name
        self.long_name = long_name
        self.units     = units
        self.mode      = mode
        self.col_id    = col_id
        self.unit_mult = unit_mult
        self.unit_off  = unit_off

    def alloc_raw(self,nraw):

        self.rawdata = np.zeros((nraw))


class timetype:

    def __init__(self,ntimes):
    
        self.year  = -9*np.ones((ntimes))
        self.month = -9*np.ones((ntimes))
        # This is a floating point decimal day
        self.day   = -9.0*np.ones((ntimes))

        # This is a decimal datenumber
        self.datenum = -9.0*np.ones((ntimes))

    def setbounds(self):

        self.yeara  = int(self.year[0])
        self.montha = int(self.month[0])
        self.yearz  = int(self.year[-1])
        self.monthz = int(self.month[-1])


# =======================================================================================

def load_xml(xmlfile):
    
    import xml.etree.ElementTree as et
    constants = []
    variables = []
    
    xmlroot = et.parse(xmlfile).getroot()
    
    for elem in xmlroot.iter('constant'):
        name      = elem.find('name').text.strip()
        long_name = elem.find('long_name').text.strip()
        units     = elem.find('units').text.strip()
        mode      = elem.find('units').text.strip()
        value     = float(elem.find('value').text)
        constants.append(contype(name,long_name,units,mode,value))


    for elem in xmlroot.iter('variable'):
        name      = elem.find('name').text.strip()
        long_name = elem.find('long_name').text.strip()
        units     = elem.find('units').text.strip() 
        mode      = elem.find('mode').text.strip()
        col_id    = int(elem.find('col_id').text)
        unit_mult = float(elem.find('unit_mult').text)
        unit_off  = float(elem.find('unit_off').text)
        variables.append(vartype(name,long_name,units,mode,col_id,unit_mult,unit_off))
    
    csv_file      = xmlroot.find('csv_file').text.strip()
    time_res_sec  = float(xmlroot.find('time_resolution').text)
    n_header_row  = int(xmlroot.find('n_header_row').text)
    n_fields      = int(xmlroot.find('n_fields').text)
    nodata_flags  = xmlroot.find('time_resolution').text.strip()
    grid_out_name = xmlroot.find('grid_name_out').text.strip()
    fill_value    = xmlroot.find('fill_value_out').text.strip()
    missing_value = xmlroot.find('missing_value_out').text.strip()

    ctrl_params=ctrltype(csv_file,time_res_sec,n_header_row,n_fields,nodata_flags, \
                grid_out_name,fill_value,missing_value)


    return(constants,variables,ctrl_params)

# =======================================================================================

def usage():
     print('')
     print('=======================================================================')
     print('')
     print(' python ConvertMetCSVtoCLMALM.py -h --f=<xml-file-name>')
     print('')
     print('')
     print(' -h --help ')
     print('     print this help message')
     print('')
     print('')
     print(' --f=<xml-file-name>')
     print('     This is the full path to the XML file that describes the conversion')
     print('     The default packaged with this script is convert_controls.xml')
     print('')
     print('')
     print('=======================================================================')


def interp_args(argv):

    argv.pop(0)  # The script itself is the first argument, forget it

    # Name of the conversion file
    xmlfile = "none"

    try:
        opts, args = getopt.getopt(argv, 'h',["f="])

    except getopt.GetoptError as err:
        print('Argument error, see usage')
        usage()
        sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o in ("--f"):
            xmlfile = a
        else:
            assert False, "unhandled option"


    if (xmlfile == "none"):
        print("You must specify an xml file, see usage:")
        usage()
        sys.exit(2)


    return (xmlfile)


# =======================================================================================

def load_csv(ctrlp,variables):

    nlines = 0
    with open(ctrlp.csv_file, 'rb') as csvfile:
        
        csvfile.seek(ctrlp.n_header_row)
        csvreader = csv.reader(csvfile, delimiter=',')
        for row in csvreader:
            nlines+=1

        nlines-=ctrlp.n_header_row

        # Allocate raw data
        for var in variables:
            var.alloc_raw(nlines)

        # Allocate timing data
        rawtime = timetype(nlines)

        # Load raw data
        csvfile.seek(0)
        csvreader = csv.reader(csvfile, delimiter=',')
        
        iidx=0
        for idx,rowtext in enumerate(csvreader):
            if (idx>=ctrlp.n_header_row):

                # Transfer CSV data into raw vectors
                for var in variables:
                    var.rawdata[iidx] = float(rowtext[var.col_id-1])

                # Timing information
                date_str1 = rowtext[0]
                date_str2 = rowtext[1]
                date1,time1 = date_str1.split(' ')
                date2,time2 = date_str2.split(' ')
                mo1,dy1,yr1 = date1.split('/')
                mo2,dy2,yr2 = date2.split('/')
                hr1,mn1 = time1.split(':')
                hr2,mn2 = time2.split(':')
                sod1 = int(hr1)*3600 + int(mn1)*60
                sod2 = int(hr2)*3600 + int(mn2)*60

                t1 = date2num(datetime(int(yr1),int(mo1),int(dy1),int(hr1),int(mn1)))
                t2 = date2num(datetime(int(yr2),int(mo2),int(dy2),int(hr2),int(mn2)))
                teff = np.mean([t1,t2])
                datestamp = num2date(teff)

                if(datestamp.month==2 & datestamp.day==29):
                    print('LEAP DAY')
                    exit(2)
                
                rawtime.datenum[iidx] = teff
                rawtime.year[iidx]  = datestamp.year
                rawtime.month[iidx] = datestamp.month
                rawtime.day[iidx]   = float(datestamp.day) + \
                              float(datestamp.hour)/24.0 + \
                              float(datestamp.minute)/1440.0 + \
                              float(datestamp.second)/86400.0
                iidx+=1

    return(variables,rawtime)

# ========================================================================================
# ========================================================================================
#                                        Main
# ========================================================================================
# ========================================================================================

def main(argv):

    # Interpret the arguments to the script
    xmlfile = interp_args(argv)
    
    constants,variables,ctrl_params = load_xml(xmlfile)


    # Algorithm:

    # 1) Loop through file and load data
    # 2) Determine set of complete months
    # 3) Loop through set of complete months
    # 4) Read in lines, parse data, convert where necessary

    # 1)
    variables,timing = load_csv(ctrl_params,variables)

    # 2) Set max/min timing info
    timing.setbounds()

    # 3) Convert datatypes
    for var in variables:
        var.rawdata = var.rawdata*var.unit_mult + var.unit_off

    # 4) Loop output files and write

    for iyr in range(timing.yeara,timing.yearz+1):
        if(iyr == timing.yeara):
            ima = timing.montha
        else:
            ima = 1
        if(iyr == timing.yearz):
            imz = timing.monthz
        else:
            imz = 12

        for imo in range(ima,imz+1):

            ncfilename = ctrl_params.grid_name_out+"."+'{:04d}'.format(iyr)+"-"+'{:02d}'.format(imo)+".nc"
            print('Preparing: '+ncfilename)

            if(imo==12):
                imo_end = 1
                iyr_end = iyr+1
            else:
                imo_end = imo+1
                iyr_end = iyr

            datenum_a = date2num(datetime(int(iyr),int(imo),int(1)))
            datenum_b = date2num(datetime(int(iyr_end),int(imo_end),int(1)))

            # Find all the time-points betwixt
            ids = np.where((timing.datenum>=datenum_a) & (timing.datenum<datenum_b))[0]

            # Open the netcdf file
            fp = netcdf.netcdf_file(ncfilename, 'w')

            fp.acknowledgement = 'Acknowledgments: This data was quality checked, cleaned and assured by Boris Faybishenko (bafaybishenko@lbl.gov). Contact for questions about formatting can be addressed to Ryan Knox (rgknox@lbl.gov)'

            fp.createDimension('time',len(ids))
            fp.createDimension('lon',1)
            fp.createDimension('lat',1)
            fp.createDimension('scalar',1)
            
            time_out    = fp.createVariable('time','f',('time',))
            time_out[:] = timing.day[ids]-1.0
            time_out.units = 'days since {:04d}'.format(iyr)+'-{:02d}-01 00:00:00'.format(imo)
            time_out.calendar = 'noleap'

            for var in variables:
                var_out = fp.createVariable(var.name,'f',('time','lat','lon'))
                var_out[:,0,0] = var.rawdata[ids]
                var_out.units = var.units
                var_out.long_name = var.long_name
                var_out.mode = var.mode
                fp.flush()


            fp.close()



#    code.interact(local=locals())





    exit(0)




# =======================================================================================
# This is the actual call to main
   
if __name__ == "__main__":
    main(sys.argv)








