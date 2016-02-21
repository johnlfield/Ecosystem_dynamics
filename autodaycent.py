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


def DDcentEVI(ddc_fpath, sch_file, target_path, id, out_files, ddclist_fpath="", site_file_out=False):
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
    if site_file_out:
        subprocess.call("%s -s %s -n %s -W %s.100" % (ddc_fpath, sch_file, id, id), shell=True)
    else:
        subprocess.call("%s -s %s -n %s" % (ddc_fpath, sch_file, id), shell=True)
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


def bioenergy_ghg_avoid_gCO2eq_m2(biomass_gC_m2_array, conversion_tech):
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
    biomass_c_conc = 0.45
    biomass_array = (biomass_gC_m2_array / biomass_c_conc) * 0.001   # kg biomass / m2 / y
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
    print "The following site / equilibrium / scenario combinations will be simulated: "
    for run_combo in run_combos:
        print "   %s / %s / %s" % (run_combo[0], run_combo[1], run_combo[2])
    proceed = raw_input("Select 'q' to quit or 'enter' to proceed: ")
    if proceed == 'q':
        return
    print

    # clear out current_results directory
    for file in glob.glob(os.path.join(recent_path, '*.*')):
        os.remove(file)

    # execute simulations
    start = time.time()
    for run_combo in run_combos:

        # read run specifications
        site = run_combo[1]
        bin_file = run_combo[2].split('.')[0]+'.bin'
        sch_file = run_combo[3]
        bin_handle = site+'_'+bin_file.split('.')[0]
        handle = bin_handle+'_'+sch_file.split('.')[0]
        print "Simulating %s / %s / %s..." % (site, bin_file, sch_file)
        print

        # define all input file paths
        schedule_fpath = base_path+sch_file
        site_path = sites_path+site+'/'
        site_fpath = site_path+bin_handle+'.100'
        weather_fpath = site_path+'weather.wth'
        soil_fpath = site_path+'soils.in'
        binary_fpath = site_path+site+'_'+bin_file
        if not os.path.exists(binary_fpath):
            print "*** ERROR: No spin-up exists at ", binary_fpath
            print
            return

        # move files, execute simulations, and clean up
        copy_in(work_path, schedule_fpath, site_fpath, weather_fpath, soil_fpath, library_path)
        os.rename(work_path+bin_handle+'.100', work_path+'site.100')
        DDcentEVI(ddc_fpath, sch_file, work_path, handle, [], ddclist_fpath=ddclist_fpath)
        clean_out(work_path, recent_path, ['.bin', '.lis'])
        print
        print

    # report out
    print "****************************************************************************"
    print
    time_sec = round((time.time() - start), 2)
    time_min = round(time_sec/60.0, 2)
    print "It took %.2f minutes total to execute the specified DayCent simulation ensemble." % time_min
    print
    return


def spinup(base_path, sites_path, work_path, library_path, ddc_fpath, ddclist_fpath, run_combos):
    """
    """
    # determine spin-ups to run, and associated input files
    spinup_combos = []
    for run_combo in run_combos:
        spinup_combos.append((run_combo[1], run_combo[2]))
    unique_spinup_combos = list(set(spinup_combos))
    unique_spinup_combos.sort()
    print "Please specify a spin-up scenario to run, or 'a' for all:"
    for i, combo in enumerate(unique_spinup_combos):
        print "   ", i, combo[0], "-", combo[1]
    index = raw_input("")
    runs = []
    if index == 'a':
        runs = unique_spinup_combos
    else:
        index = int(index)
        runs.append(unique_spinup_combos[index])

    # execute all specified spin-ups
    for run in runs:
        site = run[0]
        sch_file = run[1]
        print "Running spin-up %s for site %s..." % (sch_file, site)

        # determine all necessary paths
        site_path = sites_path+site+'/'
        site_fpath = site_path+'site.100'
        weather_fpath = site_path+'weather.wth'
        soil_fpath = site_path+'soils.in'
        schedule_fpath = base_path+sch_file

        # execute spin-up and plot som to verify equilibrium
        start = time.time()
        copy_in(work_path, schedule_fpath, site_fpath, weather_fpath, soil_fpath, library_path)
        handle = site+'_'+sch_file.split('.')[0]
        DDcentEVI(ddc_fpath, sch_file, work_path, handle, [], ddclist_fpath=ddclist_fpath, site_file_out=True)

        results = read_full_out(work_path+handle+'.lis', 3, 1)
        years = []
        for year in range(len(results[1])):
            years.append(year+1)
        plt.plot(years, list(results[1]))
        plt.xlabel('simulation time (years)')
        plt.ylabel('total SOM (gC/m2)')
        plt.savefig(handle+'.png')
        plt.close()

        # clean up and report out
        clean_out(work_path, site_path, ['.bin', '.lis', '.png', '%s.100' % handle])
        time_sec = round((time.time() - start), 2)
        time_min = round(time_sec/60.0, 2)
        print "It took %.2f minutes total to run the spin-up." % time_min
        print
        print
    return


def analysis(recent_path, archive_path, run_combos, labels):
    # create results archive
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H.%M")
    print "Please enter a short descriptive title for this analysis run, using underscores"
    description = raw_input("in place of spaces:  ")
    results_path = archive_path+timestamp+"_"+description+"/"
    os.mkdir(results_path)
    out_fpath = results_path+"Table.csv"
    header = ["Simulation", "Avg. harvestable NPP (MgC/ha/y)", "Avg. SOC change (MgC/ha/y)",
              "Avg. system C change (MgCO2eq/ha/y)", "Avg. biomass harvest (Mg biomass/ha/y)",
              "Avg. bioenergy GHGE avoidance (MgCO2eq/ha/y)", "Avg. total GHGE mitigation (MgCO2eq/ha/y)"]
    out_file_object = open(out_fpath, "wb")
    c = csv.writer(out_file_object)
    c.writerow(header)

    # initialize Fig. 3a-b
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2, sharex=ax1)
    ax3 = fig.add_subplot(2, 2, 3, sharex=ax1)

    # define structure for storing running curve minima and maxima for filling
    array_length = 72
    fill = {
            "m": [[0, 0], [0, 0], [0, 0]],
            "g": [[0, 0], [0, 0], [0, 0]],
            "r": [[0, 0], [0, 0], [0, 0]],
            "b": [[0, 0], [0, 0], [0, 0]],
            # "c": [[0, 0], [0, 0], [0, 0]]
           }
    panel12_fill = {
                    "m": True,
                    "g": True,
                    "r": False,
                    "b": True
                    }
    for key in fill:
        for panel in fill[key]:
            panel[0] = [999999] * array_length
            panel[1] = [-999999] * array_length

    # step through run_combos, read associated output file, and add traces as appropriate to Fig. 3a-b
    for run_combo in run_combos:
        system, site, eq_file, sch_file, tech, panel12, panel3, color, linestyle = run_combo
        handle = site+'_'+eq_file.split('.')[0]+'_'+sch_file.split('.')[0]
        lis_fpath = recent_path+handle+'.lis'
        summary = [handle]

        # read all results from file
        results = read_full_out(lis_fpath, 3, 1)
        years = []
        for year in range(len(results[0])):
            years.append(year+1)
        time = np.array(results[0])
        somsc = np.array(results[1])
        crmvst = np.array(results[2])
        aglivc = np.array(results[3])
        stdedc = np.array(results[4])
        fbrchc = np.array(results[5])
        rlwodc = np.array(results[6])
        wood1c = np.array(results[7])
        wood2c = np.array(results[8])
        totsysc = np.array(results[9])
        agcacc = np.array(results[10])
        fbracc = np.array(results[11])
        rlwacc = np.array(results[12])
        tcrem = np.array(results[13])

        # compute cumulative harvestable (aboveground) C production, and record average
        if system == 'grass':
            cumulative_harvestable_gC_m2 = np.cumsum(agcacc)
            cumulative_harvestable_MgC_ha = cumulative_harvestable_gC_m2 * 0.01
        elif system == 'forest':
            cumulative_harvestable_gC_m2 = fbrchc
            cumulative_harvestable_gC_m2 += rlwodc
            cumulative_harvestable_gC_m2 += wood1c
            cumulative_harvestable_gC_m2 += wood2c
            cumulative_harvestable_MgC_ha = cumulative_harvestable_gC_m2 * 0.01
        elif system == 'bm':
            cumulative_harvestable_gC_m2 = np.cumsum(crmvst)
            cumulative_harvestable_MgC_ha = cumulative_harvestable_gC_m2 * 0.01
        summary.append(round(cumulative_harvestable_MgC_ha[-1]/float(len(cumulative_harvestable_MgC_ha)), 3))

        # compute SOC and total ecosystem C, and record averages
        tot_SOC_MgC_ha = somsc * 0.01
        summary.append(round(np.mean(np.diff(tot_SOC_MgC_ha)), 4))
        tot_sysC_MgC_ha = totsysc * 0.01
        summary.append(round(np.mean(np.diff(tot_sysC_MgC_ha)) * 3.67, 3))

        # compute and plot mitigation from ecosystem C changes and bioenergy production
        initial_sysC_gC_m2 = totsysc[0]
        sysC_diff_gC_m2 = totsysc - initial_sysC_gC_m2
        sysC_diff_gC02eq_m2 = sysC_diff_gC_m2 * 3.67
        ghge_gCO2eq_m2 = sysC_diff_gC02eq_m2
        wood_harvest = np.array([0])
        wood_harvest = np.append(wood_harvest, np.diff(tcrem))
        logistic_efficiency = 0.9
        harvest_gc_m2 = (crmvst + wood_harvest) * logistic_efficiency
        nonzero_harvest_gc_m2 = []
        for harvest in harvest_gc_m2:
            if harvest > 0:
                nonzero_harvest_gc_m2.append(harvest)
        avg_harvest_gC_m2 = np.sum(nonzero_harvest_gc_m2) / float(len(nonzero_harvest_gc_m2))
        summary.append(round((avg_harvest_gC_m2*0.01)/0.45, 3))
        if tech:
            ghg_avoid_gCO2eq_m2 = bioenergy_ghg_avoid_gCO2eq_m2(harvest_gc_m2, tech)
            nonzero_ghg_avoid_gCO2eq_m2 = []
            for ghge in ghg_avoid_gCO2eq_m2:
                if ghge > 0:
                    nonzero_ghg_avoid_gCO2eq_m2.append(ghge)
            summary.append(round((np.sum(nonzero_ghg_avoid_gCO2eq_m2)/float(len(nonzero_ghg_avoid_gCO2eq_m2)))*0.01, 3))
            cum_ghg_avoid_gCO2eq_m2 = np.cumsum(ghg_avoid_gCO2eq_m2)
            ghge_gCO2eq_m2 += cum_ghg_avoid_gCO2eq_m2
        else:
            summary.append(0)
        ghge_MgCO2eq_ha = ghge_gCO2eq_m2 * 0.01
        summary.append(round(np.mean(np.diff(ghge_MgCO2eq_ha)), 3))

        # record results to file, and add to plot where appropriate
        c.writerow(summary)
        if panel12:
            ax1.plot(years, cumulative_harvestable_MgC_ha, color=color, linestyle=linestyle)
            fill[color][0][0] = np.minimum(fill[color][0][0], cumulative_harvestable_MgC_ha)
            fill[color][0][1] = np.maximum(fill[color][0][1], cumulative_harvestable_MgC_ha)
            ax2.plot(years, tot_sysC_MgC_ha, color=color, linestyle=linestyle)
            fill[color][1][0] = np.minimum(fill[color][1][0], tot_sysC_MgC_ha)
            fill[color][1][1] = np.maximum(fill[color][1][1], tot_sysC_MgC_ha)
        if panel3:
            ax3.plot(years, ghge_MgCO2eq_ha, color=color, linestyle=linestyle)
            fill[color][2][0] = np.minimum(fill[color][2][0], ghge_MgCO2eq_ha)
            fill[color][2][1] = np.maximum(fill[color][2][1], ghge_MgCO2eq_ha)

    # add fill
    for color in fill:
        if panel12_fill[color]:
            ax1.fill_between(years, fill[color][0][0], fill[color][0][1], facecolor=color, alpha=0.15, linewidth=0.0)
            ax2.fill_between(years, fill[color][1][0], fill[color][1][1], facecolor=color, alpha=0.15, linewidth=0.0)
        ax3.fill_between(years, fill[color][2][0], fill[color][2][1], facecolor=color, alpha=0.15, linewidth=0.0)

    # update plot formatting and add axis labels
    matplotlib.rcParams.update({'font.size': 10})
    ax1.set_title("Cumulative harvestable NPP", size=11, fontweight='bold')
    ax1.set_ylabel("MgC / ha")
    ax2.set_title("Current total ecosystem carbon", size=11, fontweight='bold')
    ax2.set_ylabel("MgC / ha")
    ax3.set_title("Cumulative GHGE benefit", size=11, fontweight='bold')
    ax3.set_ylabel("MgCO2eq / ha")
    ax3.set_xlabel("Time after harvest/LUC (years)", size=10)

    # add legend
    dummy_point = (-1, -1)
    for entry in labels:
        label, color, linestyle = entry
        ax3.plot(dummy_point[0], dummy_point[1], label=label, color=color, linestyle=linestyle)
    ax3.legend(bbox_to_anchor=(1.05, 1), loc=2, prop={'size': 10})

    # save figure, table, and report out
    out_file_object.close()
    shutil.copy(results_path+'Table.csv', archive_path+'Table.csv')
    plt.savefig(results_path+'Fig3.png')
    shutil.copy(results_path+'Fig3.png', archive_path+'Fig3.png')
    plt.close()
    results_path = results_path.split('/')[-2]
    print "Resulting figure and tabular output archived in the directory "+results_path
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

# define model runs in tuple format
#             (system,   site, spin-up.bin,                simulation.sch,            tech,       ax1/2,  ax3, color, linestyle)
run_combos = [
              ('forest', 'IA', 'pine_medSOCcrop_eq.sch',          'forest.sch',        '',  True, True, 'm',  '-'),  # Unmanaged forest regrowth, Iowa, young stand
              ('forest', 'IA', 'pine_medSOCcrop_forest_eq.sch',   'forest.sch',        '',  True, True, 'm', '-.'),  # Unmanaged forest regrowth, Iowa, 70yo stand
              ('forest', 'LA', 'pine_medSOCcrop_eq.sch',          'forest.sch',        '',  True, True, 'm', '--'),  # Unmanaged forest regrowth, Louisiana, young stand
              ('forest', 'LA', 'pine_medSOCcrop_forest_eq.sch',   'forest.sch',        '',  True, True, 'm',  ':'),  # Unmanaged forest regrowth, Louisiana, 70yo stand

              ('grass',  'IA', 'grass_lowSOCcrop_eq.sch',          'grass.sch',        '',  True, True, 'g',  '-'),  # Unmanaged grassland regrowth, Iowa, low-SOC soil
              ('grass',  'IA', 'grass_medSOCcrop_eq.sch',          'grass.sch',        '',  True, True, 'g', '-.'),  # Unmanaged grassland regrowth, Iowa, medium-SOC soil
              ('grass',  'LA', 'grass_lowSOCcrop_eq.sch',          'grass.sch',        '',  True, True, 'g', '--'),  # Unmanaged grassland regrowth, Louisiana, low-SOC soil
              ('grass',  'LA', 'grass_medSOCcrop_eq.sch',          'grass.sch',        '',  True, True, 'g',  ':'),  # Unmanaged grassland regrowth, Louisiana, medium-SOC soil

              ('forest', 'IA', 'pine_eq.sch',                'forest_harv.sch', 'current', False, True, 'r', '--'),  # Forest bioenergy, Iowa, current technology
              ('forest', 'LA', 'pine_eq.sch',                'forest_harv.sch', 'current', False, True, 'r',  '-'),  # Forest bioenergy, Louisiana, current technology

              ('bm',     'IA', 'grass_lowSOCcrop_eq.sch', 'switchgrass_IA.sch', 'current',  True, True, 'b',  '-'),  # Grass bioenergy, Iowa, low-SOC soil, current technology
              ('bm',     'IA', 'grass_medSOCcrop_eq.sch', 'switchgrass_IA.sch', 'current',  True, True, 'b', '-.'),  # Grass bioenergy, Iowa, medium-SOC soil, current technology
              ('bm',     'LA', 'grass_lowSOCcrop_eq.sch', 'switchgrass_LA.sch', 'current',  True, True, 'b', '--'),  # Grass bioenergy, Louisiana, low-SOC soil, current technology
              ('bm',     'LA', 'grass_medSOCcrop_eq.sch', 'switchgrass_LA.sch', 'current',  True, True, 'b',  ':'),  # Grass bioenergy, Louisiana, medium-SOC soil, current technology

              # ('grass',  'IA', 'grass_lowSOCcrop_eq.sch',         'switchgrass_IA.sch', 'mature',  False, True, 'c',  '-'),  # Grass bioenergy, Iowa, low-SOC soil, mature technology
              # ('grass',  'IA', 'grass_medSOCcrop_eq.sch',         'switchgrass_IA.sch', 'mature',  False, True, 'c', '-.'),  # Grass bioenergy, Iowa, medium-SOC soil, mature technology
              # ('grass',  'LA', 'grass_lowSOCcrop_eq.sch',         'switchgrass_LA.sch', 'mature',  False, True, 'c', '--'),  # Grass bioenergy, Louisiana, low-SOC soil, mature technology
              # ('grass',  'LA', 'grass_medSOCcrop_eq.sch',         'switchgrass_LA.sch', 'mature',  False, True, 'c',  ':')   # Grass bioenergy, Louisiana, medium-SOC soil, mature technology
              ]

labels = [('Line colors:', 'None', 'None'),
          (' Unmanaged forest regrowth', 'm', '-'),
          (' Unmanaged grassland regrowth', 'g', '-'),
          (' Forest bioenergy', 'r', '-'),
          (' Grass bioenergy, current technology', 'b', '-'),
          #('Grass bioenergy, future technology', 'c', '-'),

          ('Line styles:', 'None', 'None'),
          (' Iowa, young stand or low-SOC soil', 'k', '-'),
          (' Iowa, 70yo stand or medium-SOC soil', 'k', '-.'),
          (' Louisiana, young stand or low-SOC soil', 'k', '--'),
          (' Louisiana, 70yo stand or medium-SOC soil', 'k', ':')]

# initiate interactive routine
prompt = """
Please select from the following options:
   e = execute DayCent simulation ensemble and analyze results
   s = run DayCent model spin-up (individually or batch)
   a = analyze current simulation results set only
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
        execute(base_path, sites_path, work_path, library_path, recent_path, ddc_fpath,
                ddclist_fpath, run_combos)
        print
        print
        print "   Analyzing current simulation results set..."
        analysis(recent_path, archive_path, run_combos, labels)
        print
        print
    elif command == 's':
        print "   Running DayCent model spin-up..."
        spinup(base_path, sites_path, work_path, library_path, ddc_fpath, ddclist_fpath, run_combos)
        print
        print
    elif command == 'a':
        print "   Analyzing current simulation results set..."
        analysis(recent_path, archive_path, run_combos, labels)
        print
        print
    elif command == 'q':
        print "   Quitting..."
        print
        print
        exit()



#ToDo: update this to generate output commensurate with manuscript Table 1
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

#ToDo: this is orphaned code associated with generation bioenergy & BAU NPP, Rh, and total system C
# for site in sites:
#         print "****************************************************************************"
#         print
#         print "Conducting simulations for site '%s'..." % site
#         site_path = sites_path+site+'/'
#         site_fpath = site_path+'site.100'
#         weather_fpath = site_path+'weather.wth'
#         soil_fpath = site_path+'soils.in'
#
#         for run_combo in run_combos:
#             scenario = run_combo[0]
#             if site.startswith(scenario):
#                 # determine binary off which to extend
#                 binary_fpath = site_path+site+'_'+run_combo[1]
#                 if os.path.exists(binary_fpath):
#                     print "*** Simulating the %s scenario extending off of binary %s..." % (scenario, binary_fpath)
#                     print
#                 else:
#                     print "*** ERROR: No spin-up exists at ", binary_fpath
#                     print
#                     return
#                 # initialize dynamics detail plot
#                 f, panels = plt.subplots(3, sharex=True)
#
#                 simulations = ['BAU', 'bioenergy']
#                 colors = ['b', 'r']
#                 for i, simulation in enumerate(simulations):
#                     print "******* Executing %s simulation..." % simulation
#                     print
#                     index = i+2
#                     schedule_file = run_combo[index]
#                     schedule_fpath = base_path+schedule_file
#                     handle = site+'_'+scenario+'_'+simulation
#
#                     # move files, execute simulations, add to dynamics detail plot, and clean up
#                     copy_in(work_path, schedule_fpath, site_fpath, weather_fpath, soil_fpath, library_path)
#                     DDcentEVI(ddc_fpath, schedule_file, work_path, handle, out_files, ddclist_fpath=ddclist_fpath)
#                     print
#                     print
#                     results = read_full_out(work_path+handle+'_dc_sip.csv', 1, 0, delim='c')
#                     days = []
#                     for day in range(len(results[0])):
#                         days.append(day+1)
#                     for j, trace in enumerate([['Total C (gC/m2)', 56], ['NPP (gC/m2/d)', 31], ['Rh (gC/m2/d)', 21]]):
#                         panels[j].plot(days, results[trace[1]], label=simulation, color=colors[i])
#                         panels[j].set_ylabel(trace[0])
#                     clean_out(work_path, recent_path, ['.bin', '.lis'])
#
#                 # finalize plot
#                 plt.xlabel('simulation time (days)')
#                 plt.suptitle("\n Scenario '%s', bioenergy (%s) vs. BAU (%s) using spin-up %s" %
#                              (scenario, run_combo[3], run_combo[2], run_combo[1]))
#                 panels[0].legend(loc=4, prop={'size': 11})
#                 panels[0].set_ylim([0, 16000])
#                 panels[1].set_ylim([0, 10])
#                 panels[2].set_ylim([0, 5])
#                 # panels[2].set_xlim([0, 2000])
#                 fig = matplotlib.pyplot.gcf()
#                 fig.set_size_inches(12, 7)
#                 f.savefig(recent_path+run_combo[0]+'_dynamics_detail.png')