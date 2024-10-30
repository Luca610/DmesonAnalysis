import pandas as pd
import sys
import argparse
import yaml
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
import re
from scipy.interpolate import interp1d




def calculate_efficiencies_with_unc(df, selToApply):
    eff = len(df.query(selToApply)) / len(df)
    return eff, np.sqrt(eff * (1 - eff) / len(df))

def efficiency_plot(axis, axis2, pt_range_string, bkgBDTlims, eff_prompt_arr, s_eff_prompt_arr, eff_FD_arr, s_eff_FD_arr, proposed_cut):
    # BKG
    axis.set_title(f"BDT selection efficiencies Pt: {pt_range_string}", fontsize=18)
    axis.set_xlabel('BDT score (a.u)', fontsize=18)
    axis.set_ylabel(f'efficiency', fontsize=18)
    if proposed_cut > 0.001:
        interp = interp1d(bkgBDTlims, 0.1*eff_FD_arr[:,0] + 0.9*eff_prompt_arr[:,0], kind='cubic')
        axis.axvline(x=proposed_cut, color='b', linestyle='--', label=f'proposed cut efficiency: {interp(proposed_cut):.3f}')

    axis.set_ylim(0., 1.1)
    # Actual plotting
    axis.errorbar(bkgBDTlims, eff_prompt_arr[:,0], yerr=s_eff_prompt_arr[:,0], fmt='o', label='BDT efficiency prompt', markersize=4, color='red')
    axis.errorbar(bkgBDTlims, eff_FD_arr[:,0], yerr=s_eff_FD_arr[:,0], fmt='o', label='BDT efficiency FD', markersize=4, color='blue')
    axis.errorbar(bkgBDTlims, 0.1*eff_FD_arr[:,0] + 0.9*eff_prompt_arr[:,0] , yerr=s_eff_FD_arr[:,0], fmt='o', label='average eff', markersize=4, color='black') 
    axis.legend(fontsize=15)

    # Prompt
    axis2.set_title(f"BDT selection efficiencies Pt: {pt_range_string}", fontsize=18)
    axis2.set_xlabel('BDT score (a.u)', fontsize=18)
    axis2.set_ylabel(f'efficiency', fontsize=18)
    axis2.set_ylim(0., 1.1)
    # Actual plotting
    axis2.errorbar(bkgBDTlims, eff_prompt_arr[-1,:], yerr=s_eff_prompt_arr[-1,:], fmt='o', label='BDT efficiency prompt', markersize=4, color='red')
    axis2.errorbar(bkgBDTlims, eff_FD_arr[-1,:], yerr=s_eff_FD_arr[-1,:], fmt='o', label='BDT efficiency FD', markersize=4, color='blue')
    axis2.errorbar(bkgBDTlims, 0.1*eff_FD_arr[-1,:] + 0.9*eff_prompt_arr[-1,:] , yerr=s_eff_FD_arr[-1,:], fmt='o', label='average eff', markersize=4, color='black') 
    axis2.legend(fontsize=15)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
    parser.add_argument('--output_directory','-d',default='.',
                    help='Directory for output plots. If ".", save in the same directory of the input file.')
    args = parser.parse_args()

    # Opening configfile and sanity check
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    path_list = inputCfg['modelapplieddirs']
    file_name = inputCfg['outfile']
    cutVars = inputCfg['cutvars']
    proposed_cuts = inputCfg['proposed_cuts']
    assert len(cutVars['bkgBDT']['min']) == len(cutVars['bkgBDT']['max']) == len(path_list), f"a directory for each pt bin has not been provided"

    # Creates output folder for plots
    if args.output_directory == '.':
        root_plot_dir = os.path.join(Path(file_name), "efficiency_plots")
    else: 
        root_plot_dir = os.path.join(Path(args.output_directory), "efficiency_plots")
    if not os.path.exists(root_plot_dir):
        os.makedirs(root_plot_dir)
    
    # Initialise efficiency plot
    fig, ax = plt.subplots(nrows = 3, ncols = 4, figsize=(42,24))
    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, hspace=0.3, wspace=0.3)
    fig2, ax2 = plt.subplots(nrows = 3, ncols = 4, figsize=(42,24))
    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, hspace=0.3, wspace=0.3)
    
    pattern = r'[\d\.]+_[\d\.]+'
    for iPt, path in enumerate(path_list):
        row = int(iPt / 4)
        column = iPt % 4
        eff_prompt_arr = np.zeros([inputCfg['nsteps'],inputCfg['nsteps']])
        s_eff_prompt_arr = np.zeros([inputCfg['nsteps'],inputCfg['nsteps']])
        eff_FD_arr = np.zeros([inputCfg['nsteps'],inputCfg['nsteps']])
        s_eff_FD_arr = np.zeros([inputCfg['nsteps'],inputCfg['nsteps']])
        pt_range_string = re.search(pattern, path).group(0)
        print(f'Processing pt range {pt_range_string}')
        modelappl_prompt_path = path + 'Prompt_pT_' + pt_range_string + '_ModelApplied.parquet.gzip'
        modelappl_FD_path = path + 'FD_pT_' + pt_range_string + '_ModelApplied.parquet.gzip'
        df_prompt = pd.read_parquet(modelappl_prompt_path, engine='pyarrow')
        df_FD = pd.read_parquet(modelappl_FD_path, engine='pyarrow')
        bkgBDTlims = np.linspace(cutVars['bkgBDT']['min'][iPt], cutVars['bkgBDT']['max'][iPt], inputCfg['nsteps'] + 1)
        promptBDTlims = np.linspace(cutVars['promptBDT']['min'][iPt], cutVars['promptBDT']['max'][iPt], inputCfg['nsteps'] + 1)
        proposed_cut = proposed_cuts[iPt]
        if (inputCfg["bkgOnly"]):
            print(f"fixing prompt cut to {cutVars['promptBDT']['min'][iPt]}") 
            promptBDTlims =  promptBDTlims[0:2]
        for iBDTbkg, bkgBDTmax in enumerate(bkgBDTlims[1:]):
            for iBDTprompt, promptBDTmin in enumerate(promptBDTlims[:-1]):
                sel = f'ML_output_Bkg < {bkgBDTmax} & ML_output_Prompt > {promptBDTmin}'
                eff_prompt_arr[iBDTbkg, iBDTprompt], s_eff_prompt_arr[iBDTbkg, iBDTprompt] = calculate_efficiencies_with_unc(df_prompt, sel)
                eff_FD_arr[iBDTbkg, iBDTprompt], s_eff_FD_arr[iBDTbkg, iBDTprompt] = calculate_efficiencies_with_unc(df_FD, sel)
            # print(f'pT: {pt_range_string}, BDt Score < {bkgBDTmax} \n      prompt_eff = {eff_prompt} pm {s_eff_prompt} \n FD_eff = {eff_FD} pm {s_eff_FD} ')
        efficiency_plot(ax[row][column], ax2[row][column], pt_range_string, bkgBDTlims[1:], eff_prompt_arr, s_eff_prompt_arr, eff_FD_arr, s_eff_FD_arr, proposed_cut)
    print(f"saving {root_plot_dir}/BDT_efficiency_bkg.png")
    fig.savefig(f'{root_plot_dir}/BDT_efficiency_bkg.png')
    fig2.savefig(f'{root_plot_dir}/BDT_efficiency_prompt.png')