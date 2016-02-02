#!/bin/python
   
"""This module automates a set of DayCent model runs as specified in a user-selected
csv-formatted runtable, includes archiving and analysis of important output.  It requires
the following directory structure:

  Main directory: csv-formatted runtable
  Subdirectories:
    'scripts'- contains this script and supporting modules (test002.py, detail008.py)
    'library100'- contains archived sets of all necessary supporting DayCent .100 files, 
      switch.sch (formatted for automated modification), outfiles.in, and outvars.txt
    'results'- contains archived (timestamped) .bin, .lis, and .out output, as well as
      derivative plots and a log file
    'override'- contains archived sets of override .sch & site.100 files
    'runtable'- contains sites/treatments/yields.csv to which to compare results
    'workspace'- an empty directory into which temporary working files are copied
     
  * Note that .wth, soils.in, and site.100 files will be downloaded from their
    respective network repositories automatically; UNIX-executable DDcentEVI, DDClist100,
    and DDCsitarchive are called remotely.  The associate paths are hard-coded into this 
    script, but are logged it in the archive logfile.  
"""


import csv
import datetime
from ddc_tools import param_100
import glob
import os
import shutil
import subprocess
import sys
import time


mlras = ('1', '2', '3', '4A', '4B', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
         '20', '21', '22A', '22B', '23', '24', '25', '26', '27', '28A', '28B', '29', '30', '31', '32', '34A', '34B',
         '35', '36', '38', '39', '40', '41', '42', '43A', '43B', '43C', '44A', '44B', '46', '47', '48A', '48B', '49',
         '51', '52', '53A', '53B', '53C', '54', '55A', '55B', '55C', '56', '57', '58A', '58B', '58C', '58D', '60A',
         '60B', '61', '62', '63A', '63B', '64', '65', '66', '67A', '67B', '69', '70A', '70B', '70C', '70D', '71', '72',
         '73', '74', '75', '76', '77A', '77B', '77C', '77D', '77E', '78A', '78B', '78C', '79', '80A', '80B', '81A',
         '81B', '81C', '81D', '82A', '82B', '83A', '83B', '83C', '83D', '83E', '84A', '84B', '84C', '85', '86A', '86B',
         '87A', '87B', '88', '89', '90A', '90B', '91A', '91B', '92', '93A', '93B', '94A', '94B', '94C', '94D', '95A',
         '95B', '96', '97', '98', '99', '101', '102A', '102B', '102C', '103', '104', '105', '106', '107A', '107B',
         '108A', '108B', '108C', '108D', '109', '110', '111A', '111B', '111C', '111D', '111E', '112', '113', '114A',
         '114B', '115A', '115B', '115C', '116A', '116B', '116C', '117', '118A', '118B', '119', '120A', '120B', '120C',
         '121', '122', '123', '124', '125', '126', '127', '128', '129', '130A', '130B', '131A', '131B', '131C', '131D',
         '133A', '133B', '134', '135A', '135B', '136', '137', '138', '139', '140', '141', '142', '143', '144A', '144B',
         '145', '146', '147', '148', '149A', '149B', '150A', '150B', '151', '152A', '152B', '153A', '153B', '153C',
         '153D', '154', '155', '156A', '156B', '157', '158', '159A', '159B', '160', '161A', '161B', '162', '163', '164',
         '165', '166', '167', '190', '191', '192', '193', '194', '195', '196', '197', '220', '221', '222', '223', '224',
         '225', '226', '227', '228', '229', '230', '231', '232', '233', '234', '235', '236', '237', '238', '239', '240',
         '241', '242', '243', '244', '245', '246', '270', '271', '272', '273')


def outfiles(outfiles_fpath, lis_file):
    """Returns a list of DayCent output files being produced as a basis for concatenation of results.
    :param outfiles_fpath: file/path for outfiles.in used in DayCent runs (str)
    :return: list of unique output files to concatenate (list of str)
    """
    # read outfiles.in and make a list of all specified output files
    output_files = [lis_file]
    outfile = open(outfiles_fpath, "r")
    # skip the header line
    next(outfile)
    for line in outfile:
        if int(line.split()[0]) == 1:
            output_files.append(line.split()[1])
    outfile.close()
    return output_files


def daily_check(run_fpath):
    daily_treatments = False
    my_file = open(run_fpath, 'rU')
    run_file = csv.reader(my_file)
    columns = run_file.next()
    index = columns.index("Daily_out")
    next(run_file, None)
    for line in run_file:
        if int(line[index]) == 1:
            daily_treatments = True
    my_file.close()
    return daily_treatments


def out_reader(out_fpath, names_row, head_skip, tail_skip, study, site, treatment, cat_fpath, round,
               difference_columns=[]):
    """Reads a DayCent output file line-by-line, adds in a leading label column, and either saves it to the specified
    cat_fpath (including the associated header) or appends it to an already-existing file there
    Args-
        out_fpath (str):  file/path to output file to be read
        names_row (int):  row containing column names (counting Pythonically from zero!!)
        head_skip (int):  number of header lines to skip over before starting to read data
        tail_skip (int):  number of lines to skip over at end of file
        study, site, treatment (str):  identifiers to be included in the first entry of each row of data
        cat_fpath (str):  file/path to concatenate data to
        difference_columns (lis of str):  data columns for which to compute annual difference (added as a separate
            column at the end of the table)
    """
    out_file = open(out_fpath, 'r')
    contents = [[]]
    for line in out_file:
        data = line.split()
        data.insert(0, study)
        data.insert(1, site)
        data.insert(2, treatment)
        contents.append(data)
    del contents[0]
    out_file.close()

    # clean up headers and tailing line skips
    names = contents[names_row]
    names[0] = "Study"  # replace specific run id names with general categorical labels
    names[1] = "Site"
    names[2] = "Treatment"
    #ToDo:  note that this INT has to go back to TEXT for the calval routine to work!!!!
    types = ["TEXT", "INT", "TEXT", "REAL"]
    for k in range(len(names)-4):
        types.append("REAL")
    for i in range(head_skip):
        del contents[0]
    for j in range(tail_skip):
        del contents[-1]

    # apply float() to all non-id entries to convert any scientific notation instances to a decimal format
    for row in range(len(contents)):
        for entry in range(3, len(contents[row])):
            contents[row][entry] = float(contents[row][entry])

    # convert float date to int (rounded down)
    if round:
        for l in range(len(contents)):
            contents[l][3] = int(float(contents[l][3]))

    # compute all specified annual differences
    for column in difference_columns:
        if column in names:
            names.append("d_"+str(column))
            types.append("REAL")
            index = names.index(column)
            previous_value = contents[0][index]
            for i in range(len(contents)):
                current_value = contents[i][index]
                contents[i].append(current_value - previous_value)
                previous_value = current_value

    # append data to existing concatenated output file, or create one
    if os.path.exists(cat_fpath):
        write_file = open(cat_fpath, "a")
    else:
        write_file = open(cat_fpath, "w")
        contents.insert(0, names)
        contents.insert(1, types)
    for row in contents:
        row_string = ""
        for entry in row:
            row_string += str(entry)+","
        row_string = row_string[:-1]
        write_file.write(row_string+"\n")
    write_file.close()


def tsv_concatenate(out_files, read_path, write_path, study, site, treatment, round):
    """Concatenates tab-separated-value-formatted DayCent output files (.lis, various .out).
    Args-
        out_files (list of str):  names of all DayCent output files to be concatenated
        read_path (str):  path from which individual DayCent output files are read
        write_path (str):  path to which final concatenated output will be saved
        study, site, treatment (str):  identifiers to be included in the first entry of each row of data
        round (bool):  specifies whether or not to round years
    """
    for out_file in out_files:
        if os.path.exists(read_path+out_file):
            if out_file.endswith(".lis"):
                out_reader(read_path+out_file, 1, 4, 1, study, site, treatment, write_path+out_file, round, ["somsc"])
            else:
                out_reader(read_path+out_file, 0, 1, 0, study, site, treatment, write_path+out_file, round)


def cfarm_archive(ddcarch_fpath, cfarm_run, archive_fpath, target_path, model_version, n_lim=''):
    """Extracts a CFARM archive site.100 file to a target directory
    Args-
        ddcarch_fpath (str): file/path to UNIX-executable CFARM archive extractor program
        cfarm_run (int): CFARM archive run id #
        archive_fpath (str): path/file for the actual site100 archive .dsa file
        target_path (str): path where to extract site file (do not include 'site.100' filename)
    """
    target_fpath = target_path+"site.100"
    site_command = "%s -r %s %s %s" % (ddcarch_fpath, cfarm_run, archive_fpath, target_fpath)
    subprocess.call("%s" % site_command, shell=True)

    # routine for overwriting crazy site.100 soil mineral N values
    if n_lim:
        # open reader and writer objects
        site_reader = open(target_fpath)
        site_writer = open(target_fpath+'.temp', "wb")
        # iterate through the lines
        for i, line in enumerate(site_reader):
            if 252 < i < 263:
                entries = line.split()
                if float(entries[0]) > n_lim:
                    print "  ** NOTE: overwriting excessive site.100 file initial soil mineral N levels **"
                    entries[0] = n_lim
                    line = ' %s          %s\n' % (str(entries[0]), entries[1])
                site_writer.write(line)
            else:
                site_writer.write(line)
        # close files and rename
        site_reader.close()
        site_writer.close()
        os.remove(target_fpath)
        os.rename(target_fpath+".temp", target_fpath)

    if model_version == 'C':
        skip = range(63, 88) + range(93, 98) + range(134, 162) + range(304, 629)
        in_file = open(target_path+'site.100')
        out_file = open(target_path+'site.100.temp', 'w')
        line_count = 0
        for row in in_file:
            line_count += 1
            if line_count not in skip:
                out_file.write(row)
        in_file.close()
        out_file.close()
        os.remove(target_path+'site.100')
        os.rename(target_path+'site.100.temp', target_path+'site.100')


def soil_wth(soil_path, mukey, wth_path, narr, target_path, override_path, start_year):
    """Copies the appropriate soils.in and weather.wth files to a target directory
    Args-
        soil_path (str): path to soils.in file repository
        mukey (int): SSURGO map unit key
        wth_path (str): path to .wth file repository
        narr (str): NARR .wth file name (format 'NARR_narrx_narry.wth')
        target_path (str): path where to extract site file (do not include a file name)
        start_year- year at which the file is rearranged to begin, with earlier data
            appended to the end of the file (int)
    """
    # check for an override soils.in file
    over_soil_fpath = override_path+str(mukey)+".in"
    if os.path.exists(over_soil_fpath):
        shutil.copy(over_soil_fpath, target_path+"soils.in")
        print "Using manually-generated soil file from override archive."
    else:
        # extra logic added to handle soils.in archive in 3-digit sub-directories or lumped together in a single big
        # temporary repository
        soil_location = soil_path+mukey+".in"
        if not os.path.exists(soil_location):
            soil_folder = mukey[:-3]+"/"
            soil_location = soil_path+soil_folder+mukey+".in"
        shutil.copy(soil_location, target_path+"soils.in")

    wth_target_fpath = target_path+"weather.wth"
    shutil.copy(wth_path+narr+".wth", target_path+"weather.wth")
    # read .wth file and split into two blocks at the beginning of start_year
    # only perform these operations if start year is 2008 or earlier; otherwise final .wth file will be tiny/empty
    if start_year <= 2008:
        wth = csv.reader(open(wth_target_fpath, 'rU'), delimiter="\t")
        keep = []
        for row in wth:
            if int(row[2]) >= start_year:
                keep.append(row)
        with open(wth_target_fpath+".temp", 'w') as f_out:
            f_out.writelines('\t'.join(i) + '\n' for i in keep)
        os.remove(wth_target_fpath)
        os.rename(wth_target_fpath+".temp", wth_target_fpath)


def sch_build(ecotype, est_year, plant_doy, harv_skip, FRST, sgn1_doy, sgn2_doy, sene_doy, harv_doy, sim_years, target_path,
              sch_file):
    """Updates the wildcards in a standard sg0/1/2.sch template file within a target 
    directory as per a runtable entry and re-saves as sg.sch.  Note that the DayCent .sch 
    header block starting year must be indexed back by 1. 
    Args-
        est_year (int): year of crop establishment
        plant_doy (int): crop planting day
        harv_skip (int): number of harvests skipped during establishment phase
        sgn1_doy (int): post-establishment N application #1 day
        sgn2_doy (int): post-establishment N application #2 day
        sene_doy (int): senescence dat (in case of GDD mode this is earliest possible date)
        harv_doy (int): harvest event day
        sim_years (int): number of years to simulate
        target_path (str): path where to find sg0/1/2.sch template files and create final sg.sch
        sch_file (str): the name that the resulting schedule file will be given (incl. extension)
    """
    # define values for .sch template wildcards based on runtable contents
    osy1 = int(est_year)-1
    eot = int(est_year)+sim_years
    e_frst = int(plant_doy)+1
    cult = int(plant_doy)-1
    sgc = int(harv_doy)+1  # adding in delGro's CULT CSG event for faster deadfall
    lst = min(sgc+1, 365)
    harv_start = osy1+int(harv_skip)
    harv_end = int(harv_start)+sim_years
    crop_entry = ''
    if ecotype == 'Lowland':
        crop_entry = 'SG3L'
    elif ecotype == 'Upland':
        crop_entry = 'SG3U'

    # different template .sch files depending on how many harvests are skipped
    if int(harv_skip) == 0:
        sch_temp = 'sg0.sch'
        end = eot
    elif int(harv_skip) == 1:
        sch_temp = 'sg1.sch'
        end = eot
    else:
        sch_temp = 'sg2.sch'
        end = harv_end
    
    # open appropriate template file, do replacements, and save final schedule file
    sch_temp_fpath = target_path+sch_temp
    sch_fpath = target_path+sch_file
    infile = open(sch_temp_fpath)
    outfile = open(sch_fpath, 'w')
    replacements = {'@STY': str(est_year), '@LSY': str(end), '@EEY': str(est_year), '@ESY': str(osy1), '@PEY': str(eot),
                    '@PSY': str(est_year), '@22E': str(harv_start), '@22S': str(est_year), '@23E': str(harv_end),
                    '@23S': str(harv_start), '@PD': str(plant_doy), '@FD': str(e_frst), '@FRST':str(FRST), '@SN': str(sene_doy),
                    '@N1': str(sgn1_doy), '@N2': str(sgn2_doy), '@HV': str(harv_doy), '@CD': str(cult),
                    '@SGX': crop_entry, '@CT': str(sgc), '@LST': str(lst)}
    for row in infile:
        for src, target in replacements.iteritems():
            row = row.replace(src, target)
        outfile.write(row)
    infile.close()
    outfile.close()


def fert_build(sgn1_rate, sgn2_rate, frac_nh4, source_path, target_path):
    """Updates the wildcards in a standard fert.100 within a source directory as per a 
    runtable entry and save into a target directory.
    Args-
        sgn1_rate (float): post-establishment N application #1 rate (gN/m2)
        sgn2_rate (float): post-establishment N application #2 rate (gN/m2)
        source_path (str): path to directory from which the fert.100 template is copied
        target_path (str): path where final fert.100 file is saved
    """
    in_fpath = target_path+'fert.100'
    if os.path.exists(in_fpath):
        os.remove(in_fpath)
    shutil.copy(source_path+"fert.100", target_path+"fert.100")
    reader = open(in_fpath)
    outfile = target_path+'fert.100.temp'
    writer = open(outfile, 'w')
    frac_no3 = round(1 - float(frac_nh4), 2)
    frac_nh4 = round(float(frac_nh4), 2)
    replacements = {'@SR1': str(sgn1_rate), '@SR2': str(sgn2_rate), '@NH4': str(frac_nh4), '@NO3': str(frac_no3)}
    for row in reader:
        for src, target in replacements.iteritems():
            row = row.replace(src, target)
        writer.write(row)
    reader.close()
    writer.close()
    os.remove(in_fpath)
    os.rename(outfile, in_fpath)
    

def daycent(ddc_fpath, ddclist_fpath, sch_file, target_path, lis_name, suppress_out=False):
    """Runs the specified versions of DDcentEVI and DDClist100 for the specified schedule
    file located in the target directory, saving (renaming) all specified daily output 
    files, with all output files named as per the specified ID. Note that any .out files 
    to preserve must also be specified in the outfiles.in file in referenced library 
    archive.  Note- having trouble linking DDC programs and input files in different 
    directories unless the overall working directory is the same as that where the input 
    files are located.  
    Args-
        ddc_fpath (str):  file/path to UNIX-executable DailyDayCent program
        ddclist_fpath (str):  file/path to UNIX-executable list100 program
        sch_file (str):  filename of schedule file in target directory
        target_path (str):  directory where all model inputs are located and results are saved
        lis_name (str):  name to be given to all model output files
        suppress_out (boolean):  suppress DayCent output screen printing?  (default=False)
    """
    os.chdir(target_path)
    if suppress_out:
        null = " >/dev/null"
    else:
        null = ""
    subprocess.call("%s -s %s -n %s%s" % (ddc_fpath, sch_file, lis_name, null), shell=True)
    subprocess.call("%s %s %s outvars.txt" % (ddclist_fpath, lis_name, lis_name), shell=True)


def double_print(statement, write_file_object):
    """Simultaneously prints a statement to the screen as well as writes it as new line in the specified file object.
    Args-
        statement (str):  text to be double-printed
        write_file_object:  as defined in the calling function by "write_file_object = open(log_fpath, "w" or "a")"
    """
    print statement
    write_file_object.write(statement+"\n")


def unique(my_list):
    my_set = {}
    map(my_set.__setitem__, my_list, [])
    return my_set.keys()


def run_ddc(run_fpath,
            param_fpath="/data/paustian/AFRI/simulations/calibration/params.csv",
            param_ver='',
            mod_ver="K",
            threads=1,
            cluster="N",
            results_arch_path="/data/paustian/AFRI/simulations/results/",
            soil_path="/data/paustian/ernie/SSURGO_master_script/soil_test2/",
            wth_path="/data/paustian/AFRI/NARR_gridxy_wth/",
            site_arch_fpath="/data/paustian/AFRI/simulations/site100/runfile_cometfarm_spinup_current.dsa",
            descrip=''):
    """This function loops through a standardized runtable an executes the corresponding
    model runs.
    NOTE: parallelized version (cluster_ddc.py) only works with model version K (to reduce copy operations)
    NOTE: if you update defaults here, they have to be manually updated in cl_ddc.py as well!
    Args-
        run_fpath (str):  full file path to .csv-formatted runtable to execute
        mod_ver (str, optional):  specifies the daily DayCent model version to run, either 'K' for Ken's
            DDcentEVI (default), or 'C' for Cindy's DailyDayCent
        threads (int, optional):  total number of separate processes to run via forking (default=1)
        results_arch (str, optional):  directory to which all concatenated model output and QC
            files will be transferred upon completion (defaults to /data/paustian/...)
        cluster (Y/N, optional):  specifies special results& QC file archiving procedure needed for compatibility
            with cluster_ddc.py (default=N)
        soil_path, wth_path (str, optional):  paths to soil.in and .wth file repositories (defaults to
            /data/paustian/...)
        site_arch_fpath (str, optional):  full file path to CFARM .dsa site.100 binary archive (defaults to a
            local copy of 'runfile_cometfarm_spinup_current.dsa')
    """
    # define all /run_ddc/ sub-directory pathways, and create necessary temporary work directories
    abspath = os.path.abspath(__file__)           # get absolute path where script is located
    script_path = os.path.dirname(abspath)        # associated directory only
    os.chdir(script_path)
    os.chdir('..')
    run_ddc_path = os.getcwd()
    library_path = run_ddc_path+"/library100/"
    override_path = run_ddc_path+"/override/"
    model_path = run_ddc_path+"/model/"  # '/data/paustian/AFRI/simulations/model_versions/Yao_AET_bare/'
    # if temporary output paths already exit, terminate the program; otherwise create them
    work_path = run_ddc_path+"/workspace/"
    temp_results_path = run_ddc_path+"/temp_results/"
    temp_qc_path = temp_results_path+"QC/"
    if os.path.exists(work_path) or os.path.exists(temp_results_path) or os.path.exists(temp_qc_path):
        print
        print "Copies of the required temporary run directories already exist in /data/paustian/AFRI/simulations/run_ddc/."
        print "Please manually delete those directories and try again."
        print "Program terminating."
        print
        print
        sys.exit()
    os.makedirs(work_path)
    os.makedirs(temp_results_path)
    os.makedirs(temp_qc_path)

    # get script version
    script = os.path.basename(__file__)
    # creation = time.ctime(os.path.getctime(abspath))  Doesn't seem to work; returns same value as modification
    modification = time.ctime(os.path.getmtime(abspath))

    # update crop.100.template and like files based on specified parameter set version
    print
    version_used = param_100(param_fpath, library_path, library_path, version=param_ver)
    print

    if not descrip:
        descrip = raw_input("Please enter a run description (no spaces):  ")

    # generate an analysis time stamp, create results & QC file archive pathways
    results_path = results_arch_path
    time_stamp = datetime.datetime.now().strftime("%Y-%m-%d,%H.%M")
    descriptor = run_fpath.split("/")[-1][:-4]
    if cluster == "N":
        results_folder = time_stamp+"__"+descriptor+"__"+str(version_used)+"__"+descrip
        results_path = results_arch_path+results_folder+"/"
        if not os.path.exists(results_path):
            os.makedirs(results_path)
        qc_path = results_path+"QC/"
    elif cluster == "Y":
        results_path = results_arch_path
        qc_path = results_arch_path[:-8]+"QC/"
    if not os.path.exists(qc_path):
        os.makedirs(qc_path)

    # specification of model version and supporting programs
    ddc_fpath = ""
    ddclist_fpath = ""
    if mod_ver == 'K' or mod_ver is None:
        ddc_fpath = model_path+"DDcentEVI"
        ddclist_fpath = model_path+"DDClist100"
    elif mod_ver == 'C':
        ddc_fpath = "/data/paustian/AFRI/simulations/model_versions/cindy_DDC/DailyDayCent"
        ddclist_fpath = "/data/paustian/AFRI/simulations/model_versions/cindy_DDC/DailyDayCent_list100"
    ddcarch_fpath = model_path+"DDCsitarchive"

    # general analysis specifications
    number_processes = int(threads)
    sch_file = 'sg.sch'
    # define .out files to preserve
    out_files = outfiles(library_path+"outfiles.in", "X.lis")
    # read the runtable and check for daily_out=1
    daily_treatments = daily_check(run_fpath)
    daily_out_files = outfiles(library_path+"outfiles_daily.in", "Y.lis")
    qc_points = []
    # qc_points = [['Anderson13', 'Urbana_1', 'SG1'],
    #              ['Bonin14', 'Fostoria', 'B50'],
    #              ['Bonin14', 'WashingtonCourtHouse', 'B50'],
    #              ['Bonin14', 'Wilkesville', 'B50'],
    #              ['Dohleman12', 'SouthFarms', 'N00'],
    #              ['Frank04', 'North', 'sites'],
    #              ['Frank04', 'South', 'sites'],
    #              ['Wilson14', 'Uthe', 'W0'],
    #              ['Wilson14', 'Uthe', 'W01'],
    #              ['Nikiema11', 'MSU_Catham', 'MSU56'],
    #              ['Schmer12', 'Mandan', 'S67'],
    #              ['Hong12', 'Bristol_1', 'H56']]
    sim_years = 30      # total number of post-establishment years to simulate for each entry

    # print, log, and archive all analysis metadata
    log_file = "%s_%s_log.txt" % (time_stamp, descriptor)
    c = open(results_path+log_file, "w")
    print
    double_print("Simulating the runtable %s" % run_fpath, c)
    double_print("Running code %s last modified %s" % (script, modification), c)
    double_print("Executed from %s across %i process threads" % (run_ddc_path, number_processes), c)
    double_print("Results archived at %s" % results_path, c)
    double_print("", c)
    double_print("Simulation specifications - ", c)
    double_print("DayCent version/path: %s" % ddc_fpath, c)
    double_print("List100 utility version/path: %s" % ddclist_fpath, c)
    double_print("Site100 archive utility version/path: %s" % ddcarch_fpath, c)
    double_print("CFARM site.100 archive version/path: %s" % site_arch_fpath, c)
    double_print("Soils.in file repository: %s" % soil_path, c)
    double_print("Weather file repository: %s" % wth_path, c)
    c.close()
    print

    # start analysis run timer, open the runtable file, set up a placeholder list of lists for job sets to be
    # distributed across the different threads, and propagate it
    start = time.time()
    run_table = open(run_fpath, 'rU')
    lines = csv.reader(run_table)
    lines.next()  # skip two headers rows
    lines.next()
    row = []
    job_sets = [[]]
    treat_count = 0
    for j in range(number_processes-1):
        job_sets.append([])
    counter = 0
    for line in lines:
        treat_count += 1
        for j in range(len(line)):
            row.append(line[j])
        job_sets[counter].append(row)
        row = []
        if counter < (number_processes-1):
            counter += 1
        else:
            counter = 0
    run_table.close()
    print "Simulating a total of %i DayCent model runs." % treat_count
    print

    # spin off the desired number of child processes, note the pid, and pass an entry from the job set list of lists
    # using http://www.petercollingridge.co.uk/blog/running-multiple-processes-python as a template
    children = []
    for process in range(number_processes):
        child = process+1  # number child processes starting with 1 (not 0!)
        time.sleep(0.1)    # throttling back on the spin-offs so printed status updates don't get mixed up
        pid = os.fork()
        if pid:            # keep track of child processes in parent process
            children.append(pid)
        else:              # within each child process, execute its set of the full job list
            print
            print "Child process #%i (PID %s) to execute %i DayCent model runs." % (child, os.getpid(),
                                                                                    len(job_sets[process]))
            print

            # create child-process-specific sub-directories within workspace/ and temp_results/
            work_path_child = work_path+str(child)+"/"
            temp_results_path_child = temp_results_path+str(child)+"/"
            os.makedirs(work_path_child)
            os.makedirs(temp_results_path_child)

            # read through the job set sub-list by sub-list
            for k in range(len(job_sets[process])):
                cfarm_run, mukey, narrx, narry, lat, longitude, study, treat, ecotype, est_year, plant_doy, sgne_rate,\
                sgne_doy, harv_skip, sgn1_rate, sgn1_doy, sgn2_rate, sgn2_doy, sene_doy, har_doy, daily_out, frac_nh4, \
                response, site, state, calval, FRST, SENM = job_sets[process][k]

                # define the basic run parameters
                run_id = '%s_%s_%s' % (study, site, treat)
                narr = "NARR_%s_%s" % (str(narrx), str(narry))
                print "Running DayCent for treatment %s at %s, %s as per %s." % (treat, site, state, study)
                print "Using a soils.in file for SSURGO mukey %s, %s.wth, and CFARM run index %s." % (str(mukey), narr,
                                                                                                      cfarm_run)

                # move entire library100 archive contents to workspace directory
                for f in glob.glob(os.path.join(library_path, '*')):
                    shutil.copy(f, work_path_child)

                # retrieve site.100 file from override archive or using cfarm_archive()
                site_fpath = override_path+site+".100"
                if os.path.exists(site_fpath):
                    shutil.copy(site_fpath, work_path_child+"site.100")
                    print "Using manually-generated site.100 file from override archive."
                else:
                    cfarm_archive(ddcarch_fpath, cfarm_run, site_arch_fpath, work_path_child, mod_ver, n_lim=20.0)
                    print "Using automatically-generated site.100 file from the CFARM archive."

                # retrieve schedule file from override archive or using sch_build()
                over_sch_fpath = override_path+treat+".sch"
                if os.path.exists(over_sch_fpath):
                    shutil.copy(over_sch_fpath, work_path_child+sch_file)
                    print "Using manually-generated schedule file from override archive."
                else:
                    sch_build(ecotype, est_year, plant_doy, harv_skip, FRST, sgn1_doy, sgn2_doy, SENM, har_doy,
                              sim_years, work_path_child, sch_file)
                    print "Using automatically-generated schedule file based on runtable."
                # raw_input("Paused...")
                print

                # retrieve soil.in & .wth, update fert.100, copy over standard outfiles.in, and run DayCent
                soil_wth(soil_path, mukey, wth_path, narr, work_path_child, override_path, int(est_year))
                fert_build(sgn1_rate, sgn2_rate, frac_nh4, library_path, work_path_child)
                if not os.path.exists(work_path_child+"outfiles.in"):
                    shutil.copy(library_path+"outfiles.in", work_path_child)
                # suppress DayCent output if in multiple-process run mode
                if number_processes == 1:
                    daycent(ddc_fpath, ddclist_fpath, sch_file, work_path_child, "X")
                else:
                    daycent(ddc_fpath, ddclist_fpath, sch_file, work_path_child, "X", suppress_out=True)
                print

                # daily-output scenario additional runs
                if bool(int(daily_out)):
                    # switch the outfiles.in in the working directory
                    os.remove(work_path_child+"outfiles.in")
                    shutil.copy(library_path+"outfiles_daily.in", work_path_child+"outfiles.in")
                    # re-run DayCent with new outfiles.in, saving output with new name (suppressing output if necessary)
                    if number_processes == 1:
                        daycent(ddc_fpath, ddclist_fpath, sch_file, work_path_child, "Y")
                    else:
                        daycent(ddc_fpath, ddclist_fpath, sch_file, work_path_child, "Y", suppress_out=True)
                    print

                # raw_input("Paused...")
                # concatenate all model output into child-specific temporary results directory
                tsv_concatenate(out_files, work_path_child, temp_results_path_child, study, site, treat, True)
                if bool(int(daily_out)):
                    tsv_concatenate(daily_out_files, work_path_child, temp_results_path_child, study, site, treat, False)

                # if the run is included in the QC list, copy all files in working directory to a new QC sub-directory
                if [study, site, treat] in qc_points:
                    qc_id_path = temp_qc_path+"%s_%s_%s/" % (study, site, treat)
                    os.makedirs(qc_id_path)
                    for each_file in glob.glob(os.path.join(work_path_child, "*")):
                        shutil.copy(each_file, qc_id_path)

                # delete all run-specific input and output files from the workspace directory
                run_files = [work_path_child+"site.100", work_path_child+"soils.in", work_path_child+"weather.wth",
                             work_path_child+sch_file, work_path_child+"fert.100", work_path_child+"X.bin",
                             work_path_child+"Y.bin", work_path_child+"outfiles.in", work_path_child+"crop.100",
                             work_path_child+"fix.100", work_path_child+"harv.100"]
                for out_file in out_files:
                    run_files.append(work_path_child+out_file)
                for daily_out_file in daily_out_files:
                    run_files.append(work_path_child+daily_out_file)
                for each in run_files:
                    if os.path.exists(each):
                        os.remove(each)

            print
            print "Child #%i done running." % child
            print
            os._exit(0)

    # wait for each child process to finish execution, and report its final 'waitpid' status
    for process, pid in enumerate(children):
        child = process+1
        results = os.waitpid(pid, 0)
        print "Child %i 'waitpid' status:" % child, results

    # calculate total and average simulation times
    sec = round((time.time() - start), 2)
    sec_per_treat = round(sec/treat_count, 2)
    minutes = round(sec/60.0, 2)
    print
    print
    print "Analysis complete."
    print "It took %s minutes total to run the %s treatments (%s sec/treatment)." % (str(minutes), str(treat_count),
                                                                                     str(sec_per_treat))

    # concatenate results from each child's temporary results directory into the timestamped results archive
    # merge out file lists into a singe list with unique entries
    if daily_treatments:
        out_files += daily_out_files
    unique_out_files = unique(out_files)
    for out_file in unique_out_files:
        i = 0
        for j in range(len(children)):
            child = j+1
            temp_results_path_child = temp_results_path+str(child)+"/"
            if os.path.exists(temp_results_path_child+out_file):
                if i == 0:
                    shutil.copy(temp_results_path_child+out_file, results_path)
                    i += 1
                else:
                    input_file = open(temp_results_path_child+out_file, 'r')
                    write_file = open(results_path+out_file, "a")
                    input_file.next()
                    input_file.next()
                    for line in input_file:
                        write_file.write(line)
                    input_file.close()
                    write_file.close()

    # move QC files to results archive, remove temp_results and working directories
    # raw_input("Paused....")
    qc_folders = os.listdir(temp_qc_path)
    print qc_folders
    for folder in qc_folders:
        shutil.move(temp_qc_path+folder, qc_path)
    shutil.rmtree(temp_results_path)
    shutil.rmtree(work_path)
    shutil.move(library_path+'crop.100', results_path+'crop.100')
    shutil.move(library_path+'fix.100', results_path+'fix.100')
    shutil.move(library_path+'harv.100', results_path+'harv.100')
    print
    print
    return results_path, version_used, descrip


def soil_wth2(soil_path, mukey, wth_path, narr, target_path, override_path):
    """Copies the appropriate soils.in and weather.wth files to a target directory
    Args-
        soil_path (str): path to soils.in file repository
        mukey (int): SSURGO map unit key
        wth_path (str): path to .wth file repository
        narr (str): NARR .wth file name (format 'NARR_narrx_narry.wth')
        target_path (str): path where to extract site file (do not include a file name)
        start_year- year at which the file is rearranged to begin, with earlier data
            appended to the end of the file (int)
    """
    # check for an override soils.in file
    over_soil_fpath = override_path+str(mukey)+".in"
    if os.path.exists(over_soil_fpath):
        shutil.copy(over_soil_fpath, target_path+"soils.in")
        print "Using manually-generated soil file from override archive."
    else:
        # extra logic added to handle soils.in archive in 3-digit sub-directories or lumped together in a single big
        # temporary repository
        soil_location = soil_path+mukey+".in"
        if not os.path.exists(soil_location):
            soil_folder = mukey[:-3]+"/"
            soil_location = soil_path+soil_folder+mukey+".in"
        shutil.copy(soil_location, target_path+"soils.in")

    shutil.copy(wth_path+narr+".wth", target_path+"weather.wth")


def daycent1(ddc_fpath, sch1, target_path, lis_name, suppress_out=False):
    """Runs the specified versions of DDcentEVI and DDClist100 for the specified schedule
    file located in the target directory, saving (renaming) all specified daily output
    files, with all output files named as per the specified ID. Note that any .out files
    to preserve must also be specified in the outfiles.in file in referenced library
    archive.  Note- having trouble linking DDC programs and input files in different
    directories unless the overall working directory is the same as that where the input
    files are located.
    Args-
        ddc_fpath (str):  file/path to UNIX-executable DailyDayCent program
        ddclist_fpath (str):  file/path to UNIX-executable list100 program
        sch_file (str):  filename of schedule file in target directory
        target_path (str):  directory where all model inputs are located and results are saved
        lis_name (str):  name to be given to all model output files
        suppress_out (boolean):  suppress DayCent output screen printing?  (default=False)
    """
    os.chdir(target_path)
    if suppress_out:
        null = " >/dev/null"
    else:
        null = ""
    subprocess.call("%s -s %s -n %s%s" % (ddc_fpath, sch1, lis_name, null), shell=True)


def daycent2(ddc_fpath, ddclist_fpath, sch2, target_path, old_lis_name, new_lis_name, suppress_out=False):
    os.chdir(target_path)
    if suppress_out:
        null = " >/dev/null"
    else:
        null = ""
    subprocess.call("%s -s %s -e %s.bin -n %s%s" % (ddc_fpath, sch2, old_lis_name, new_lis_name, null), shell=True)
    subprocess.call("%s %s %s outvars.txt" % (ddclist_fpath, new_lis_name, new_lis_name), shell=True)


def sch_build2(sch_file, first_year, sgn1_doy, sene_doy, harv_doy, last_year, target_path):
    """Updates the wildcards in a standard sg0/1/2.sch template file within a target
    directory as per a runtable entry and re-saves as sg.sch.  Note that the DayCent .sch
    header block starting year must be indexed back by 1.
    Args-
        est_year (int): year of crop establishment
        plant_doy (int): crop planting day
        harv_skip (int): number of harvests skipped during establishment phase
        sgn1_doy (int): post-establishment N application #1 day
        sgn2_doy (int): post-establishment N application #2 day
        sene_doy (int): senescence dat (in case of GDD mode this is earliest possible date)
        harv_doy (int): harvest event day
        sim_years (int): number of years to simulate
        target_path (str): path where to find sg0/1/2.sch template files and create final sg.sch
        sch_file (str): the name that the resulting schedule file will be given (incl. extension)
    """
    # define values for .sch template wildcards based on runtable contents
    lst = min(int(harv_doy)+1, 365)

    # open appropriate template file, do replacements, and save final schedule file
    infile = open(target_path+sch_file+'.sch_temp')
    outfile = open(target_path+sch_file, 'w')
    replacements = {'@STY': str(first_year), '@LSY': str(last_year), ' @N1': str(sgn1_doy), '@SN': str(sene_doy),
                    '@HV': str(harv_doy), '@LST': str(lst)}
    for row in infile:
        for src, target in replacements.iteritems():
            row = row.replace(src, target)
        outfile.write(row)
    infile.close()
    outfile.close()


def run_ddc2(run_fpath,
             param_fpath="/data/paustian/AFRI/simulations/calibration/params.csv",
             param_ver='',
             mod_ver="K",
             threads=1,
             cluster="N",
             results_arch_path="/data/paustian/AFRI/simulations/results/",
             soil_path="/data/paustian/ernie/SSURGO_master_script/soil_test2/",
             wth_path="/data/paustian/AFRI/NARR_gridxy_wth/",
             site_arch_fpath="/data/paustian/AFRI/simulations/site100/runfile_cometfarm_spinup_current.dsa",
             descrip=''):
    """This function loops through a standardized runtable an executes the corresponding
    model runs.
    NOTE: parallelized version (cluster_ddc.py) only works with model version K (to reduce copy operations)
    NOTE: if you update defaults here, they have to be manually updated in cl_ddc.py as well!
    Args-
        run_fpath (str):  full file path to .csv-formatted runtable to execute
        mod_ver (str, optional):  specifies the daily DayCent model version to run, either 'K' for Ken's
            DDcentEVI (default), or 'C' for Cindy's DailyDayCent
        threads (int, optional):  total number of separate processes to run via forking (default=1)
        results_arch (str, optional):  directory to which all concatenated model output and QC
            files will be transferred upon completion (defaults to /data/paustian/...)
        cluster (Y/N, optional):  specifies special results& QC file archiving procedure needed for compatibility
            with cluster_ddc.py (default=N)
        soil_path, wth_path (str, optional):  paths to soil.in and .wth file repositories (defaults to
            /data/paustian/...)
        site_arch_fpath (str, optional):  full file path to CFARM .dsa site.100 binary archive (defaults to a
            local copy of 'runfile_cometfarm_spinup_current.dsa')
    """
    # define all /run_ddc/ sub-directory pathways, and create necessary temporary work directories
    abspath = os.path.abspath(__file__)           # get absolute path where script is located
    script_path = os.path.dirname(abspath)        # associated directory only
    os.chdir(script_path)
    os.chdir('..')
    run_ddc_path = os.getcwd()
    library_path = run_ddc_path+"/library100/"
    override_path = run_ddc_path+"/override/"
    model_path = run_ddc_path+"/model/"  # '/data/paustian/AFRI/simulations/model_versions/Yao_AET_bare/'
    # if temporary output paths already exit, terminate the program; otherwise create them
    work_path = run_ddc_path+"/workspace/"
    temp_results_path = run_ddc_path+"/temp_results/"
    temp_qc_path = temp_results_path+"QC/"
    if os.path.exists(work_path) or os.path.exists(temp_results_path) or os.path.exists(temp_qc_path):
        print
        print "Copies of the required temporary run directories already exist in /data/paustian/AFRI/simulations/run_ddc/."
        print "Please manually delete those directories and try again."
        print "Program terminating."
        print
        print
        sys.exit()
    os.makedirs(work_path)
    os.makedirs(temp_results_path)
    os.makedirs(temp_qc_path)

    # get script version
    script = os.path.basename(__file__)
    # creation = time.ctime(os.path.getctime(abspath))  Doesn't seem to work; returns same value as modification
    modification = time.ctime(os.path.getmtime(abspath))

    # update crop.100.template and like files based on specified parameter set version
    print
    version_used = param_100(param_fpath, library_path, library_path, version=param_ver)
    print

    if not descrip:
        descrip = raw_input("Please enter a run description (no spaces):  ")

    # generate an analysis time stamp, create results & QC file archive pathways
    results_path = results_arch_path
    time_stamp = datetime.datetime.now().strftime("%Y-%m-%d,%H.%M")
    descriptor = run_fpath.split("/")[-1][:-4]
    if cluster == "N":
        results_folder = time_stamp+"__"+descriptor+"__"+str(version_used)+"__"+descrip
        results_path = results_arch_path+results_folder+"/"
        if not os.path.exists(results_path):
            os.makedirs(results_path)
        qc_path = results_path+"QC/"
    elif cluster == "Y":
        results_path = results_arch_path
        qc_path = results_arch_path[:-8]+"QC/"
    if not os.path.exists(qc_path):
        os.makedirs(qc_path)

    # specification of model version and supporting programs
    ddc_fpath = ""
    ddclist_fpath = ""
    if mod_ver == 'K' or mod_ver is None:
        ddc_fpath = model_path+"DDcentEVI"
        ddclist_fpath = model_path+"DDClist100"
    elif mod_ver == 'C':
        ddc_fpath = "/data/paustian/AFRI/simulations/model_versions/cindy_DDC/DailyDayCent"
        ddclist_fpath = "/data/paustian/AFRI/simulations/model_versions/cindy_DDC/DailyDayCent_list100"
    ddcarch_fpath = model_path+"DDCsitarchive"

    # general analysis specifications
    number_processes = int(threads)
    # define .out files to preserve
    out_files = outfiles(library_path+"outfiles.in", "Y.lis")
    qc_points = []
    # qc_points = [['Anderson13', 'Urbana_1', 'SG1'],
    #              ['Bonin14', 'Fostoria', 'B50'],
    #              ['Bonin14', 'WashingtonCourtHouse', 'B50'],
    #              ['Bonin14', 'Wilkesville', 'B50'],
    #              ['Dohleman12', 'SouthFarms', 'N00'],
    #              ['Frank04', 'North', 'sites'],
    #              ['Frank04', 'South', 'sites'],
    #              ['Wilson14', 'Uthe', 'W0'],
    #              ['Wilson14', 'Uthe', 'W01'],
    #              ['Nikiema11', 'MSU_Catham', 'MSU56'],
    #              ['Schmer12', 'Mandan', 'S67'],
    #              ['Hong12', 'Bristol_1', 'H56']]
    start_year = 1979
    conversion_year = 2015
    end_year = 2045

    # print, log, and archive all analysis metadata
    log_file = "%s_%s_log.txt" % (time_stamp, descriptor)
    c = open(results_path+log_file, "w")
    print
    double_print("Simulating the runtable %s" % run_fpath, c)
    double_print("Running code %s last modified %s" % (script, modification), c)
    double_print("Executed from %s across %i process threads" % (run_ddc_path, number_processes), c)
    double_print("Results archived at %s" % results_path, c)
    double_print("", c)
    double_print("Simulation specifications - ", c)
    double_print("DayCent version/path: %s" % ddc_fpath, c)
    double_print("List100 utility version/path: %s" % ddclist_fpath, c)
    double_print("Site100 archive utility version/path: %s" % ddcarch_fpath, c)
    double_print("CFARM site.100 archive version/path: %s" % site_arch_fpath, c)
    double_print("Soils.in file repository: %s" % soil_path, c)
    double_print("Weather file repository: %s" % wth_path, c)
    c.close()
    print

    # start analysis run timer, open the runtable file, set up a placeholder list of lists for job sets to be
    # distributed across the different threads, and propagate it
    start = time.time()
    run_table = open(run_fpath, 'rU')
    lines = csv.reader(run_table)
    lines.next()  # skip two headers rows
    lines.next()
    row = []
    job_sets = [[]]
    treat_count = 0
    for j in range(number_processes-1):
        job_sets.append([])
    counter = 0
    for line in lines:
        treat_count += 1
        for j in range(len(line)):
            row.append(line[j])
        job_sets[counter].append(row)
        row = []
        if counter < (number_processes-1):
            counter += 1
        else:
            counter = 0
    run_table.close()
    print "Simulating a total of %i DayCent model runs." % treat_count
    print

    # spin off the desired number of child processes, note the pid, and pass an entry from the job set list of lists
    # using http://www.petercollingridge.co.uk/blog/running-multiple-processes-python as a template
    children = []
    for process in range(number_processes):
        child = process+1  # number child processes starting with 1 (not 0!)
        time.sleep(0.1)    # throttling back on the spin-offs so printed status updates don't get mixed up
        pid = os.fork()
        if pid:            # keep track of child processes in parent process
            children.append(pid)
        else:              # within each child process, execute its set of the full job list
            print
            print "Child process #%i (PID %s) to execute %i DayCent model runs." % (child, os.getpid(),
                                                                                    len(job_sets[process]))
            print

            # create child-process-specific sub-directories within workspace/ and temp_results/
            work_path_child = work_path+str(child)+"/"
            temp_results_path_child = temp_results_path+str(child)+"/"
            os.makedirs(work_path_child)
            os.makedirs(temp_results_path_child)

            # read through the job set sub-list by sub-list
            for k in range(len(job_sets[process])):
                cfarm_run, mukey, narrx, narry, irr, nlcd, crp, study, treat, sgn1_rate, sgn1_doy, sene_doy, harv_doy, \
                frac_nh4, site, state, fips, sch1, sch2 = job_sets[process][k]

                # define the basic run parameters
                run_id = '%s_%s_%s' % (study, site, treat)
                narr = "NARR_%s_%s" % (str(narrx), str(narry))
                print "Running DayCent for treatment %s at %s, %s as per %s." % (treat, site, state, study)
                print "Using a soils.in file for SSURGO mukey %s, %s.wth, and CFARM run index %s." % (str(mukey), narr,
                                                                                                      cfarm_run)

                # move entire library100 archive contents to workspace directory
                for f in glob.glob(os.path.join(library_path, '*')):
                    shutil.copy(f, work_path_child)

                # retrieve site.100 file from override archive or using cfarm_archive()
                cfarm_archive(ddcarch_fpath, cfarm_run, site_arch_fpath, work_path_child, mod_ver, n_lim=20.0)

                # retrieve soil.in & .wth, update fert.100, copy over standard outfiles.in, and run DayCent
                soil_wth2(soil_path, mukey, wth_path, narr, work_path_child, override_path)
                fert_build(sgn1_rate, 0, frac_nh4, library_path, work_path_child)
                if not os.path.exists(work_path_child+"outfiles.in"):
                    shutil.copy(library_path+"outfiles.in", work_path_child)

                # retrieve schedule file from override archive or using sch_build()
                sch_build2(sch1, start_year-1, sgn1_doy, sene_doy, harv_doy, conversion_year-2, work_path_child)
                # raw_input("Paused...")

                # suppress DayCent output if in multiple-process run mode
                if number_processes == 1:
                    daycent1(ddc_fpath, sch1, work_path_child, "X")
                    os.remove(work_path_child+sch1)
                    sch_build2(sch2, conversion_year, sgn1_doy, sene_doy, harv_doy, end_year, work_path_child)
                    # raw_input("Paused...")
                    daycent2(ddc_fpath, ddclist_fpath, sch2, work_path_child, "X", "Y")
                else:
                    daycent1(ddc_fpath, sch1, work_path_child, "X", suppress_out=True)
                    os.remove(work_path_child+sch1)
                    sch_build2(sch2, conversion_year, sgn1_doy, sene_doy, harv_doy, end_year, work_path_child)
                    # raw_input("Paused...")
                    daycent2(ddc_fpath, ddclist_fpath, sch2, work_path_child, "X", "Y", suppress_out=True)
                print

                # raw_input("Paused...")
                # concatenate all model output into child-specific temporary results directory
                tsv_concatenate(out_files, work_path_child, temp_results_path_child, study, site, treat, True)

                # if the run is included in the QC list, copy all files in working directory to a new QC sub-directory
                if [study, site, treat] in qc_points:
                    qc_id_path = temp_qc_path+"%s_%s_%s/" % (study, site, treat)
                    os.makedirs(qc_id_path)
                    for each_file in glob.glob(os.path.join(work_path_child, "*")):
                        shutil.copy(each_file, qc_id_path)

                # delete all run-specific input and output files from the workspace directory
                run_files = [work_path_child+"site.100", work_path_child+"soils.in", work_path_child+"weather.wth",
                             work_path_child+sch1, work_path_child+sch2, work_path_child+"fert.100", work_path_child+"X.bin",
                             work_path_child+"Y.bin", work_path_child+"outfiles.in", work_path_child+"crop.100",
                             work_path_child+"fix.100", work_path_child+"harv.100"]
                for out_file in out_files:
                    run_files.append(work_path_child+out_file)
                for each in run_files:
                    if os.path.exists(each):
                        os.remove(each)

            print
            print "Child #%i done running." % child
            print
            os._exit(0)

    # wait for each child process to finish execution, and report its final 'waitpid' status
    for process, pid in enumerate(children):
        child = process+1
        results = os.waitpid(pid, 0)
        print "Child %i 'waitpid' status:" % child, results

    # calculate total and average simulation times
    sec = round((time.time() - start), 2)
    sec_per_treat = round(sec/treat_count, 2)
    minutes = round(sec/60.0, 2)
    print
    print
    print "Analysis complete."
    print "It took %s minutes total to run the %s treatments (%s sec/treatment)." % (str(minutes), str(treat_count),
                                                                                     str(sec_per_treat))

    # concatenate results from each child's temporary results directory into the timestamped results archive
    # merge out file lists into a singe list with unique entries
    unique_out_files = unique(out_files)
    for out_file in unique_out_files:
        i = 0
        for j in range(len(children)):
            child = j+1
            temp_results_path_child = temp_results_path+str(child)+"/"
            if os.path.exists(temp_results_path_child+out_file):
                if i == 0:
                    shutil.copy(temp_results_path_child+out_file, results_path)
                    i += 1
                else:
                    input_file = open(temp_results_path_child+out_file, 'r')
                    write_file = open(results_path+out_file, "a")
                    input_file.next()
                    input_file.next()
                    for line in input_file:
                        write_file.write(line)
                    input_file.close()
                    write_file.close()

    # move QC files to results archive, remove temp_results and working directories
    # raw_input("Paused....")
    qc_folders = os.listdir(temp_qc_path)
    print qc_folders
    for folder in qc_folders:
        shutil.move(temp_qc_path+folder, qc_path)
    shutil.rmtree(temp_results_path)
    shutil.rmtree(work_path)
    shutil.move(library_path+'crop.100', results_path+'crop.100')
    shutil.move(library_path+'fix.100', results_path+'fix.100')
    shutil.move(library_path+'harv.100', results_path+'harv.100')
    print
    print
    return results_path, version_used, descrip