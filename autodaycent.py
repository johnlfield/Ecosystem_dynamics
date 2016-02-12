#!/bin/python

"""****************************************************************************
This interactive python routine facilitates execution of DayCent simulations
and results analysis for bioenergy ecosystem dynamics assessments. It is
currently organized around the generation of Figures 3a-c for the manuscript
"Greenhouse Gas Emission Reduction from Biofuels and Ecosystem Carbon Storage",
and structured such that:
   * model execution and results analysis are de-coupled for convenience
   * ensemble runs across multiple sites are accommodated
   * analysis results are automatically archived with full provenance
To add additional sites to the model run ensemble, create a new site directory
with weather.wth and soils.in files and use the
****************************************************************************
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
import numpy as np


def gm2_to_Mgha(data_list):
    for i, datum in enumerate(data_list):
        data_list[i] = float(datum) * 0.01
    return data_list


def DDcentEVI(ddc_fpath, sch_file, target_path, id, out_files, ddclist_fpath="", bin_fpath="", site_file=False):
    """Runs the specified versions of DDcentEVI and DDClist100 for the specified schedule
    file located in the target directory, saving (renaming) all specified daily output
    files, with all output files named as per the specified ID. Note that any .out files
    to preserve must also be specified in the outfiles.in file in referenced library
    archive.
    Args:
        ddc_fpath- file/path to UNIX-executable DailyDayCent program (str)
        ddclist_fpath- file/path to UNIX-executable list100 program (str)
        sch_file- filename of schedule file in target directory (str)
        target_path- directory where all model inputs are located and results are saved (str)
        id- name to be given to all model output files (str)
        out_files- list of daily output .out files to save (list of str)
    """
    os.chdir(target_path)
    if bin_fpath:
        bin_file = bin_fpath.split('/')[-1]
        if site_file:
            subprocess.call("%s -s %s -n %s -i %s -W site.100 >/dev/null" % (ddc_fpath, sch_file, id, bin_file), shell=True)
        else:
            subprocess.call("%s -s %s -n %s -i %s >/dev/null" % (ddc_fpath, sch_file, id, bin_file), shell=True)
    else:
        if site_file:
            subprocess.call("%s -s %s -n %s -W site.100 >/dev/null" % (ddc_fpath, sch_file, id), shell=True)
        else:
            subprocess.call("%s -s %s -n %s >/dev/null" % (ddc_fpath, sch_file, id), shell=True)
    if ddclist_fpath:
        subprocess.call("%s %s %s outvars.txt" % (ddclist_fpath, id, id), shell=True)
    for file in out_files:
        os.rename(target_path+file, target_path+id+"_"+file)


def auto_type(datum):
    """Convert data types to integer or float if possible
    """
    try:
        return int(datum)
    except:
        try:
            return float(datum)
        except:
            return str(datum)


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


def read_full_out(out_fpath, head_skip, tail_skip, delim='s'):
    """Reads space-delimited .lis or .out output files (including skipping headers and
    dummy head or tail data rows), and returns averages of entire columns or averages of
    annual differences within a column.
    Args:
        out_fpath- full path and filename of output file to be analyzed (str)
        head_skip- number of rows to skip at beginning of file, incl. headers and dummy data (int)
        tail_skip- numbers of rows to skip at end of file, typically dummy data (int)
    """
    # read all lines into a 2D list, delete the list initialization and specified head/tail lines
    input = open(out_fpath, 'rU')
    mylist = [[]]
    if delim == 's':
        for i in range(head_skip):
            next(input)
        for line in input:
            mylist.append(line.split())
        del mylist[0]
        for j in range(tail_skip):
            del mylist[-1]

    else:
        if delim == 't':
            lines = csv.reader(input, delimiter="\t")
        elif delim == 'c':
            lines = csv.reader(input)
        for line in lines:
            mylist.append(line)
        del mylist[0]
        for i in range(head_skip):
            del mylist[0]
        for j in range(tail_skip):
            del mylist[-1]

    # convert all entries in the 2D list to float format where possible, or zero in the case
    # of very small numbers recorded in scientific notation
    for k in range(len(mylist)):
        for l in range(len(mylist[k])):
            mylist[k][l] = auto_type(mylist[k][l])

    #return transpose
    mylist = zip(*mylist)
    return mylist


def yie_ghg_sum(id, work_path):
    """For the specified bioenergy or BAU scenario, this function takes basic output file
    column-averaged output from read_out(), and performs the appropriate unit conversions
    to report scenario biomass yield, grain yield, change in SOM, N2O emissions, and total
    CO2 equivalent GHGs.
    Args:
        id- string describing simulating, incl. soil, LUC scenario, bioenergy/BAU (str)
    """
    #call read_out() to analyze the simulation .lis and year_summary.out output
    avg_cols_out1, diff_cols_out1 = summarize_out(work_path+id+'.lis', 4, 1, [6,8], [7,9])
    avg_cols_out2, diff_cols_out2 = summarize_out(work_path+id+'_year_summary.out', 1, 0, [1], [])
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


def bioenergy_ghg_avoid_gCO2eq_m2(biomass_gC_m2_list, conversion_tech):
    """Based on JohnF_BiomassBioenergy_2_9_16.xlsx
    """
    # define conversion technology-dependant constants
    if conversion_tech == 'current':
        biofuel_yield = 9.441   # MJ/kg biomass, cell C24
        biopower_yield = 1.274   # MJ/kg biomass, cell C25
        ghg_int_electricity = 166.6   # gCO2eq/MJ, cell D55
    elif conversion_tech == 'mature':
        biofuel_yield = 12.194   # MJ/kg biomass, cell C36
        biopower_yield = 0.230   # MJ/kg biomass, cell C37
        ghg_int_electricity = 142.6   # gCO2eq/MJ, cell D55

    # compile GHG avoidance factors that reflect lifecycle emissions of fossil fuel baseline and bioenergy supply chain
    ghg_int_biofuel_supplychain = 9.43   # gCO2eq/MJ, cell C20
    ghg_int_gasoline = 92.40   # gCO2eq/MJ, cell D54
    biofuel_ghg_avoid_MJ = ghg_int_gasoline - ghg_int_biofuel_supplychain   # gCO2eq/MJ
    biopower_ghg_avoid_MJ = ghg_int_electricity - ghg_int_biofuel_supplychain   # gCO2eq/MJ

    # calculate energy production and GHG avoidance from biomass supply array
    biomass_c_array = np.array(biomass_gC_m2_list)   # gC / m2 / y
    biomass_c_conc = 0.45
    biomass_array = (biomass_c_array / biomass_c_conc) * 0.001   # kg biomass / m2 / y
    biofuel_ghg_avoid_array = biomass_array * biofuel_yield * biofuel_ghg_avoid_MJ   # gCO2eq / m2 / y
    biopower_ghg_avoid_array = biomass_array * biopower_yield * biopower_ghg_avoid_MJ   # gCO2eq / m2 / y
    bioenergy_ghg_avoid_array = biofuel_ghg_avoid_array + biopower_ghg_avoid_array   # gCO2eq / m2 / y

    return bioenergy_ghg_avoid_array


def copy_in(work_path, schedule_fpath, site_fpath, weather_fpath, soil_fpath, library_path, bin_fpath=""):
    """Clear out workspace, then copy in schedule, weather, soil & library files
    """
    for file in glob.glob(os.path.join(work_path, '*.*')):
        os.remove(file)
    shutil.copy(site_fpath, work_path)
    shutil.copy(schedule_fpath, work_path)
    shutil.copy(weather_fpath, work_path)
    shutil.copy(soil_fpath, work_path)
    if bin_fpath:
        shutil.copy(bin_fpath, work_path)
    for file in glob.glob(os.path.join(library_path, '*')):
        shutil.copy(file, work_path)
    return


def clean_out(work_path, move_path, save_types):
    """Archive key results and clean up working directory
    """
    for file in glob.glob(os.path.join(work_path, '*')):
        move_status = False
        for save_type in save_types:
            if file.endswith(save_type):
                move_status = True
        if move_status:
            file_name = file.split('/')[-1]
            move_fpath = move_path+file_name
            shutil.move(file, move_fpath)
        else:
            os.remove(file)
    return


def execute(base_path, sites_path, work_path, library_path, recent_path, ddc_fpath, ddclist_fpath, run_combos):
    # determine all model runs and prompt user to proceed
    print "The following equilibrium/BAU/bioenergy file combinations will be executed: "
    for run_combo in run_combos:
        print "   %s:  %s / %s / %s" % (run_combo[0], run_combo[1], run_combo[2], run_combo[3])
    sites = [x.split('/')[-1] for x in os.walk(sites_path).next()[1]]
    print "Simulations will be executed across the following sites: ", sites
    print
    proceed = raw_input("Select 'q' to quit or 'enter' to proceed: ")
    if proceed == 'q':
        return
    print

    # clear out current_results directory
    for file in glob.glob(os.path.join(recent_path, '*.*')):
        os.remove(file)

    # execute simulations
    start = time.time()
    out_files = ['year_summary.out', 'dc_sip.csv']
    for site in sites:
        print "****************************************************************************"
        print
        print "Conducting simulations for site '%s'..." % site
        site_path = sites_path+site+'/'
        site_fpath = site_path+'site.100'
        weather_fpath = site_path+'weather.wth'
        soil_fpath = site_path+'soils.in'

        for run_combo in run_combos:
            scenario = run_combo[0]
            if site.startswith(scenario):
                # determine binary off which to extend
                binary_fpath = site_path+site+'_'+run_combo[1]
                if os.path.exists(binary_fpath):
                    print "*** Simulating the %s scenario extending off of binary %s..." % (scenario, binary_fpath)
                    print
                else:
                    print "*** ERROR: No spin-up exists at ", binary_fpath
                    print
                    return
                # initialize dynamics detail plot
                f, panels = plt.subplots(3, sharex=True)

                simulations = ['BAU', 'bioenergy']
                colors = ['b', 'r']
                for i, simulation in enumerate(simulations):
                    print "******* Executing %s simulation..." % simulation
                    print
                    index = i+2
                    schedule_file = run_combo[index]
                    schedule_fpath = base_path+schedule_file
                    handle = site+'_'+scenario+'_'+simulation

                    # move files, execute simulations, add to dynamics detail plot, and clean up
                    copy_in(work_path, schedule_fpath, site_fpath, weather_fpath, soil_fpath, library_path)
                    DDcentEVI(ddc_fpath, schedule_file, work_path, handle, out_files, ddclist_fpath=ddclist_fpath)
                    print
                    print
                    results = read_full_out(work_path+handle+'_dc_sip.csv', 1, 0, delim='c')
                    days = []
                    for day in range(len(results[0])):
                        days.append(day+1)
                    for j, trace in enumerate([['Total C (gC/m2)', 56], ['NPP (gC/m2/d)', 31], ['Rh (gC/m2/d)', 21]]):
                        panels[j].plot(days, results[trace[1]], label=simulation, color=colors[i])
                        panels[j].set_ylabel(trace[0])
                    clean_out(work_path, recent_path, ['.bin', '.lis'])

                # finalize plot
                plt.xlabel('simulation time (days)')
                plt.suptitle("\n Scenario '%s', bioenergy (%s) vs. BAU (%s) using spin-up %s" %
                             (scenario, run_combo[3], run_combo[2], run_combo[1]))
                panels[0].legend(loc=4, prop={'size': 11})
                panels[0].set_ylim([0, 16000])
                panels[1].set_ylim([0, 10])
                panels[2].set_ylim([0, 5])
                # panels[2].set_xlim([0, 2000])
                fig = matplotlib.pyplot.gcf()
                fig.set_size_inches(12, 7)
                f.savefig(recent_path+run_combo[0]+'_dynamics_detail.png')

    # report out
    print "****************************************************************************"
    print
    time_sec = round((time.time() - start), 2)
    time_min = round(time_sec/60.0, 2)
    print "It took %.2f minutes total to execute the specified DayCent simulation ensemble." % time_min
    print
    return


def spinup(base_path, sites_path, work_path, library_path, ddc_fpath, ddclist_fpath):
    """
    """
    # determine site for which to run spin-up, and associated input files
    sites = [x.split('/')[-1] for x in os.walk(sites_path).next()[1]]
    print "Please specify site number for which to run spin-up:"
    for i, site in enumerate(sites):
        print "   ", i, site
    index = raw_input("")
    site = sites[int(index)]
    site_path = sites_path+site+'/'
    print "Running spin-up for site '"+site+"'..."
    site_fpath = site_path+'site.100'
    weather_fpath = site_path+'weather.wth'
    soil_fpath = site_path+'soils.in'
    print "Please specify a schedule file for the spin-up:"
    schedule = raw_input("")
    schedule_fpath = base_path+schedule

    # execute spin-up and plot som to verify equilibrium
    start = time.time()
    copy_in(work_path, schedule_fpath, site_fpath, weather_fpath, soil_fpath, library_path)
    handle = site+'_'+schedule.split('.')[0]
    DDcentEVI(ddc_fpath, schedule, work_path, handle, [], ddclist_fpath=ddclist_fpath, site_file=True)

    results = read_full_out(work_path+handle+'.lis', 3, 1)
    years = []
    for year in range(len(results[0])):
        years.append(year+1)
    plt.plot(years, results[0])
    plt.xlabel('simulation time (years)')
    plt.ylabel('total SOM (gC/m2)')
    plt.savefig(handle+'.png')
    plt.close()

    # clean up and report out
    clean_out(work_path, site_path, ['.bin', '.lis', '.png', 'site.100'])
    time_sec = round((time.time() - start), 2)
    time_min = round(time_sec/60.0, 2)
    print "It took %.2f minutes total to run the spin-up." % time_min
    return


def analysis(recent_path, archive_path):
    # create results archive
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H.%M")
    print "Please enter a short descriptive title for this analysis run, using underscores"
    description = raw_input("in place of spaces:  ")
    results_path = archive_path+timestamp+"_"+description+"/"
    os.mkdir(results_path)

    # initialize Fig. 3a-b
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2, sharex=ax1)
    ax3 = fig.add_subplot(2, 2, 3, sharex=ax1)

    # define future bioenergy crop productivity factor
    future_productivity_increase = 2.0

    # search through all files in /current_results, and add traces as appropriate to Fig. 3a-b
    suffixes = ['grass_bioenergy.lis', 'grass_BAU.lis', 'forest_BAU.lis']
    labels = ['Perennial energy crops', 'Unmanaged grassland', 'Unmanaged forest']
    # label_triggers = [0, 0, 0]
    colors = ['b', 'g', 'm']
    for file in glob.glob(os.path.join(recent_path, '*')):
        for i, suffix in enumerate(suffixes):
            if file.endswith(suffix):
                file_name = file.split('/')[-1]
                results = read_full_out(file, 3, 1)
                years = []
                for year in range(len(results[0])):
                    years.append(year+1)

                # compute and plot total aboveground C
                crmvst = np.array(results[2])
                aglivc = np.array(results[3])
                stdedc = np.array(results[4])
                fbrchc = np.array(results[5])
                rlwodc = np.array(results[6])
                wood1c = np.array(results[7])
                wood2c = np.array(results[8])
                totagc = crmvst + aglivc
                # totagc += stdedc
                totagc += fbrchc
                totagc += rlwodc
                totagc += wood1c
                totagc += wood2c
                totagc.tolist()
                agcacc = results[10]
                tcrem = results[13]

                # compute and plot cumulative aboveground C production (panel 1)
                #ToDo: figure out if these two calculations are comparable...
                if suffix.startswith('grass'):
                    cumulative_harvestable = np.cumsum(agcacc).tolist()
                    cumulative_harvestable = gm2_to_Mgha(cumulative_harvestable)
                    ax1.plot(years, cumulative_harvestable, label=file_name, color=colors[i])
                elif suffix.startswith('forest'):
                    cumulative_harvestable = gm2_to_Mgha(totagc)
                    ax1.plot(years, cumulative_harvestable, label=file_name, color=colors[i])

                # compute and plot total ecosystem C (panel 2)
                totsysc = gm2_to_Mgha(list(results[9]))
                ax2.plot(years, totsysc, label=labels[i], color=colors[i])

                # compute and plot bioenergy scenario carbon benefits (panel 3)
                #ToDo: standardize colors for biomass type / LCA combos, and linestyles for different site / initial condition combos
                #ToDo: change looping structure to specifiy these within tuples- perhaps one control loop for each level
                #ToDo: add fill between extremes of site / initial condition combos
                #ToDo: get rid of the IF statements as possible; might need to specify that no calculation is done in
                # the 0 biomass yield case for computational efficiency
                if suffix.endswith('BAU.lis'):
                    initial_C = totsysc[0]
                    totsysc[:] = [((x - initial_C) * 3.67) for x in totsysc]   # gCO2eq/m2
                    ax3.plot(years, totsysc, label=file_name, color=colors[i])
                if suffix.endswith('bioenergy.lis'):
                    if suffix.startswith('grass'):
                        initial_C = totsysc[0]
                        totsysc[:] = [((x - initial_C) * 3.67) for x in totsysc]   # gCO2eq/m2
                        harvest_gc_m2 = np.array(agcacc) + np.array(tcrem)
                        linestyles = ['--', '-.']
                        for j, tech in enumerate(['current', 'mature']):
                            if tech == 'mature':
                                harvest_gc_m2 *= future_productivity_increase
                            ghg_avoid_gCO2eq_m2 = bioenergy_ghg_avoid_gCO2eq_m2(harvest_gc_m2, tech)
                            cumulative_ghg_avoid = np.cumsum(ghg_avoid_gCO2eq_m2)
                            total_ghg_mitigation = cumulative_ghg_avoid + np.array(totsysc[:])   # gCO2eq/m2
                            total_ghg_mitigation *= 0.01   # MgCO2eq/ha
                            ax3.plot(years, total_ghg_mitigation, label=file_name+'-'+tech, color=colors[i],
                                     linestyle=linestyles[j])

    # # finalize plot
    matplotlib.rcParams.update({'font.size': 10})
    ax1.set_title("Cumulative harvestable NPP", size=11, fontweight='bold')
    ax1.set_ylabel("MgC / ha")
    ax2.set_title("Current total ecosystem carbon", size=11, fontweight='bold')
    ax2.set_ylabel("MgC / ha")
    ax3.set_title("Cumulative GHGE benefit", size=11, fontweight='bold')
    ax3.set_ylabel("MgCO2eq / ha")
    ax3.set_xlabel("Time after harvest/LUC (years)", size=10)
    ax3.legend(bbox_to_anchor=(1.05, 1), loc=2, prop={'size': 10})
    plt.savefig(results_path+'Fig3.png')
    shutil.copy(results_path+'Fig3.png', archive_path+'Fig3.png')
    plt.close()
    results_path = results_path.split('/')[-2]
    print "All resulting data & metadata is archived in the directory "+results_path
    return


#########################################################################################

print
print __doc__

# define key paths
base_path = os.path.dirname(os.path.abspath(__file__))+'/'
sites_path = base_path+"sites/"
work_path = base_path+"workspace/"
library_path = base_path+"input_files/"
ddc_fpath = base_path+"model/DDcentEVI"
ddclist_fpath = base_path+"model/DDClist100"
recent_path = base_path+"current_results/"
archive_path = base_path+"results_archive/"

# define analysis scenarios
run_combos = [('forest', 'pine_eq.bin', 'pine.sch', 'pine_rot.sch'),
              ('grass', 'graze_eq.bin', 'grass.sch', 'switchgrass.sch')]

# initiate interactive routine
prompt = """
Please select from the following options:
   e = execute DayCent simulation ensemble
   s = run DayCent model spin-up
   a = analyze current simulation results set
   q = quit
"""
command = ""
while True:
    while True:
        command = raw_input(prompt)
        if command in ('e', 's', 'a', 'q'):
            break
        else:
            print "   Please try again!"
    if command == 'e':
        print "   Executing DayCent simulation ensemble..."
        print
        execute(base_path, sites_path, work_path, library_path, recent_path, ddc_fpath,
                ddclist_fpath, run_combos)
    elif command == 's':
        print "   Running DayCent model spin-up..."
        print
        spinup(base_path, sites_path, work_path, library_path, ddc_fpath, ddclist_fpath)
    elif command == 'a':
        print "   Analyzing current simulation results set..."
        print
        analysis(recent_path, archive_path)
    elif command == 'q':
        print "   Quitting..."
        print
        print
        exit()





#ToDO: update this to generate output commensurate with manuscript Table 1
# #report analysis and log analysis runtime, results summary, and results archive location
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