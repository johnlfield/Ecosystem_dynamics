#!/bin/python

"""This module automates sets of DayCent simulation runs, including basic
analysis of results and data/metadata archiving.  It will run the specified
.sch file for all combinations of soils and land use histories for a given
site for which a site.100 file exits.
    NOTE: don't rename any of the files or sub-directories or it will break!!

"""

from run_ddc import out_reader
import os
import csv
import glob
import shutil
import subprocess
import time
import datetime
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt


def DDcentEVI(ddc_fpath, ddclist_fpath, sch_file, target_path, id, out_files):
    """Runs the specified versions of DDcentEVI and DDClist100 for the specified schedule
    file located in the target directory, saving (renaming) all specified daily output
    files, with all output files named as per the specified ID. Note that any .out files
    to preserve must also be specified in the outfiles.in file in referenced library
    archive.  Note- having trouble linking DDC programs and input files in different
    directories unless the overall working directory is the same as that where the input
    files are located.
    Args:
        opsys- operating system, either 'l', 'm', or 'w' (str)
        ddc_fpath- file/path to UNIX-executable DailyDayCent program (str)
        ddclist_fpath- file/path to UNIX-executable list100 program (str)
        sch_file- filename of schedule file in target directory (str)
        target_path- directory where all model inputs are located and results are saved (str)
        id- name to be given to all model output files (str)
        out_files- list of daily output .out files to save (list of str)
    """
    os.chdir(target_path)
    subprocess.call("%s -s %s -n %s" % (ddc_fpath, sch_file, id), shell=True)
    subprocess.call("%s %s %s outvars.txt" % (ddclist_fpath, id, id), shell=True)
    for file in out_files:
        os.rename(target_path+file, target_path+id+"_"+file)


def summarize_out(out_fpath, head_skip, tail_skip, avg_cols, diff_cols):
    """Reads space-delimited .lis or .out output files (including skipping headers and
    dummy head or tail data rows), and returns averages of entire columns or averages of
    annual differences within a column.
    Args:
        out_fpath- full path and filename of output file to be analyzed (str)
        head_skip- number of rows to skip at beginning of file, incl. headers and dummy data (int)
        tail_skip- numbers of rows to skip at end of file, typically dummy data (int)
        avg_cols- columns for which to return averages, counting from 0 (list of int)
        diff_cols- columns for which to return (list of int)
    """
    #skip the proper number of headers, read each subsequent line into a 2D list, assuming
    #space-delimited, and delete the 2D list initialization row and any dummy last lines
    cont = [[]]
    input = open(out_fpath, 'r')
    for i in range(head_skip):
        next(input)
    for line in input:
        cont.append(line.split())
    del cont[0]
    for j in range(tail_skip):
        del cont[-1]

    #convert all entries in the 2D list to float format where possible, or zero in the case
    #of very small numbers recorded in scientific notation
    for k in range(len(cont)):
        for l in range(len(cont[k])):
            cont[k][l] = float(cont[k][l])

    #loop through and compute each avg_cols entry
    avg_cols_out = []
    for m in avg_cols:
        tot = 0
        count = 0
        for p in range(len(cont)):
            tot += cont[p][m]
            count += 1
        avg_cols_out.append(tot/count)

    #loop through and compute each diff_cols entry
    diff_cols_out = []
    for n in diff_cols:
        tot = 0
        count = 0
        for p in range(len(cont)-1):
            tot += cont[p+1][n]-cont[p][n]
            count += 1
        diff_cols_out.append(tot/count)

    return avg_cols_out, diff_cols_out


def read_full_out(out_fpath, head_skip, tail_skip):
    """Reads space-delimited .lis or .out output files (including skipping headers and
    dummy head or tail data rows), and returns averages of entire columns or averages of
    annual differences within a column.
    Args:
        out_fpath- full path and filename of output file to be analyzed (str)
        head_skip- number of rows to skip at beginning of file, incl. headers and dummy data (int)
        tail_skip- numbers of rows to skip at end of file, typically dummy data (int)
        avg_cols- columns for which to return averages, counting from 0 (list of int)
        diff_cols- columns for which to return (list of int)
    """
    #skip the proper number of headers, read each subsequent line into a 2D list, assuming
    #space-delimited, and delete the 2D list initialization row and any dummy last lines
    # cont = [[]]
    input = open(out_fpath, 'rU')
    lines = csv.reader(input)
    mylist = [[]]
    for line in lines:
        mylist.append(line)
    del mylist[0]
    for i in range(head_skip):
        del mylist[0]
    for j in range(tail_skip):
        del mylist[-1]

    # convert all entries in the 2D list to float or int format where possible, or zero in the case
    # of very small numbers recorded in scientific notation
    # for k in range(len(mylist)):
    #     for l in range(len(mylist[k])):
    #         for fn in (float, int):
    #             try:
    #                 mylist[k][l] = fn(mylist[k][l])
    #             except ValueError:
    #                 pass

    #return transpose
    mylist = zip(*mylist)
    return mylist


def yie_ghg_sum(id):
    """For the specified bioenergy or BAU scenario, this function takes basic output file
    column-averaged output from read_out(), and performs the appropriate unit conversions
    to report scenario biomass yield, grain yield, change in SOM, N2O emissions, and total
    CO2 equivalent GHGs.
    Args:
        id- string describing simulating, incl. soil, LUC scenario, bioenergy/BAU (str)
    """
    #call read_out() to analyze the simulation .lis and year_summary.out output
    avg_cols_out1, diff_cols_out1 = summarize_out(dirwork+id+'.lis', 4, 1, [6,8], [7,9])
    avg_cols_out2, diff_cols_out2 = summarize_out(dirwork+id+'_year_summary.out', 1, 0, [1], [])
    cmrvst = avg_cols_out1[0]   #gC/m2/y
    tcrem = diff_cols_out1[0]    #gC/m2/y
    cgrain = avg_cols_out1[1]   #gC/m2/y
    d_somsc = diff_cols_out1[1]   #gC/m2/y
    n2oflux = avg_cols_out2[0]   #gN2O-N/m2/y

    #unit conversions & data recording
    hbmyield = (cmrvst/0.45)/100.0  #Mg-BM/ha/y
    wbmyield = (tcrem/0.45)/100.0  #Mg-BM/ha/y
    gryield = (cgrain/0.45)/100.0  #Mg-grain/ha/y
    d_som = d_somsc/100.0  #Mg/ha/y
    n2o = n2oflux*(44.0/28.0)*10000.0*.001
    ghg = (-1*d_som)+(298.0*n2o*0.001)

    case = id.split('_')[-1]
    sum = [case, round(hbmyield,2), round(wbmyield,2), round(gryield,2), round(d_som,4),
           round(n2o,3), round(ghg,3)]
    return sum


#########################################################################################


#specification of schedule file name and local directory structure
print
print __doc__

# while True:
#     try:
#         site_abrv = raw_input("Please specify your site, f=florida, a=alabama, l=louisiana, i=iowa, t=test:  ")
#     finally:
#         if site_abrv in ('f', 'a', 'l', 'i', 't'):
#             break
#         else:
#             print "  Please try again!"
# if site_abrv == 'f':
#     site = 'florida'
# elif site_abrv == 'a':
#     site = 'alabama'
# elif site_abrv == 'i':
#     site = 'iowa'
# elif site_abrv == 'l':
#     site = 'louisiana'
# else:
site = 'test'
print "Please enter a short descriptive title for this analysis run, using underscores"
descrip = raw_input("in place of spaces:  ")
bau_sch_file = raw_input("Please enter the name of the BAU schedule file to run (incl. extension):  ")
bio_sch_file = raw_input("Please enter the name of the bioenergy harvest schedule file to run (incl. extension):  ")
tstamp = datetime.datetime.now().strftime("%Y-%m-%d_%H.%M")     #timestamp for archive
abspath = os.path.abspath(__file__)         #get absolute path where script is located
out_files = ['year_summary.out', 'dc_sip.csv']

script = abspath.split('/')[-1]
dirmain = os.path.dirname(abspath)+"/"      #associated directory only
dirsite = dirmain+"sites/"+site+"/"
dirwork = dirmain+"workspace/"
dirlib = dirmain+"input_files/"
dirres = dirmain+"results/"+tstamp+"_"+descrip+"/"
ddc_fpath = dirmain+"model/DDcentEVI"
ddclist_fpath = dirmain+"model/DDClist100"
print


#start analysis run timer, clear workspace, copy in schedule, .wth, & library files
start = time.time()
for file in glob.glob(os.path.join(dirwork, '*.*')):
    os.remove(file)
shutil.copy(bau_sch_file, dirwork)
shutil.copy(bio_sch_file, dirwork)
for file in glob.glob(os.path.join(dirsite, '*.wth')):
    shutil.copy(file, dirwork+"weather.wth")
for file in glob.glob(os.path.join(dirlib, '*')):
    shutil.copy(file, dirwork)


#compute soils/LUC combinations
combos = []
for s in glob.glob(os.path.join(dirsite, '*.100')):
    site_file = s.split('/')[-1]
    combos.append(site_file)
print "The following soil-LUC scenarios will be simulated for both bioenergy and BAU scenarios:  "
combos = sorted(combos)
print combos
time.sleep(3)
print
print


#for each soil/LUC combination, move the associated soils.in and site.100 file into the
#workscape, execute the simulation, convert results to .lis format, and analyze the results
stats = [[]]
for i in range(len(combos)):

    #replace the soils.in and site.100 in the workspace directory
    site_file = combos[i]
    soil_file = combos[i][:-4].split('_')[-1]+".in"
    evt_file = combos[i].split('_')[0]+"_"+combos[i].split('_')[1]+".evt"
    if os.path.exists(dirwork+"soils.in"):
        os.remove(dirwork+"soils.in")
    if os.path.exists(dirwork+"site.100"):
        os.remove(dirwork+"site.100")
    shutil.copy(dirsite+site_file, dirwork+"site.100")
    shutil.copy(dirsite+soil_file, dirwork+"soils.in")
    shutil.copy(dirsite+evt_file, dirwork)

    #set up plot
    f, panels = plt.subplots(3, sharex=True)

    #execute model runs in a loop, using the schedule-soil-luc combination as a unique identifier
    id = site_file[:-4]
    bio_id = id+"_bioenergy"
    bau_id = id+"_BAU"
    print
    for case in [[bau_sch_file, bau_id, 'g'], [bio_sch_file, bio_id, 'r']]:
        print "Running schedule file %s..." % case[0]
        DDcentEVI(ddc_fpath, ddclist_fpath, case[0], dirwork, case[1], out_files)
        # sums = yie_ghg_sum(bio_id)
        results = read_full_out(case[1]+'_dc_sip.csv', 1, 0)
        days = []
        for day in range(len(results[0])):
            days.append(day+1)
        for j, trace in enumerate([['Total C (gC/m2)', 56], ['NPP (gC/m2/d)', 31], ['Rh (gC/m2/d)', 21]]):
            panels[j].plot(days, results[trace[1]], label=case[1], color=case[2])   #, linewidth=0.6)
            panels[j].set_ylabel(trace[0])
        # sums.insert(0, combos[i].split('_')[2][:-4])
        # sums.insert(1, combos[i].split('_')[1])
        # stats.append(sums)
        print

    plt.xlabel('simulation time (days)')
    plt.suptitle('\n%s vs. %s, using spin-up %s' % (bio_sch_file, bau_sch_file, combos[i]))
    panels[0].legend(loc=4, prop={'size': 11})
    panels[0].set_ylim([0, 16000])
    panels[1].set_ylim([0, 10])
    panels[2].set_ylim([0, 5])
    # panels[2].set_xlim([0, 2000])
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(12, 7)
    f.savefig(tstamp+'_'+descrip+'.png')

# del stats[0]
print
print



# #report analysis and log analysis runtime, results summary, and results archive location
# sec = round((time.time() - start), 2)
# min = round(sec/60.0, 2)
# resfile = tstamp+"_"+descrip+"_results.csv"
# c = open(dirwork+resfile, "w")
# c.write("Analysis timstamp:  "+tstamp+'\n')
# c.write("Simulation description:  "+descrip+'\n')
# c.write("Model run automation code version:  "+script+'\n')
# c.write("DDC version:  "+ddc_fpath+'\n')
# c.write("DDClis version:  "+ddclist_fpath+'\n\n')
# c.write('soil,initial land use,scenario,herb.biomass.yield.avg (Mg/ha/y),wood.biomass.yield.avg (Mg/ha/y),\
# grain.yield.avg (Mg/ha/y),d.SOM.avg (Mg-C/ha/y),N2O.avg (kg-N2O/ha/y),GHGs.avg (Mg CO2eq/ha/y)\n')
# print
# print "*********************** ANALYSIS COMPLETE *********************** "
# print
# print "It took %s minutes total to execute and analyze these %s scenarios." % \
#        (str(min), str(2*len(combos)))
# print
# print "Logging all results & metadata as "+resfile
# print
# print "Summary of %s analysis results - " % bio_sch_file
# print
# for j in range(len(stats)):
#     print """Soil %s, '%s' initial land use, %s scenario:
#    Avg herbacious biomass yield = %s (Mg/ha/y)
#    Avg woody biomass yield = %s (Mg/ha/y)
#    Avg grain yield = %s (Mg/ha/y)
#    Avg change in SOM = %s (Mg-C/ha/y)
#    Avg N2O = %s (kg-N2O/ha/y)
#    Avg GHG footprint = %s (Mg CO2eq/ha/y)\n""" \
#              % (stats[j][0], stats[j][1], stats[j][2], stats[j][3], stats[j][4], stats[j][5],
#                 stats[j][6], stats[j][7], stats[j][8])
#     line = "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
#              % (stats[j][0], stats[j][1], stats[j][2], stats[j][3], stats[j][4], stats[j][5],
#                 stats[j][6], stats[j][7], stats[j][8])
#     c.write(line)
# c.close()

results = dirres.split('/')[-2]
print "All resulting data & metadata is archived in the directory "+results
print
print


#archive results, working directory cleanup
# shutil.move(dirwork+logfile, dirres)
if not os.path.exists(dirres):
    os.mkdir(dirres)
for file in glob.glob(os.path.join(dirwork, '*')):
    if file.endswith(".in") or file.endswith(".txt") or file.endswith(".100") or file.endswith(".wth"):
        os.remove(file)
for file in glob.glob(os.path.join(dirwork, '*')):
    if file.endswith(".bin") or file.endswith(".lis") or file.endswith(".out") or \
       file.endswith(".csv") or file.endswith(".png")or file.endswith(".sch") or file.endswith(".evt"):
        shutil.move(file, dirres)

