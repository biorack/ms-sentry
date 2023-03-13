#handle the plotting functions

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import os

pd.options.mode.chained_assignment = None  # default='warn'

def _make_unlabeled_m1_compound_figure(all_ms1_data, y_axis, compound_name, logy=False):
    ms1_data = all_ms1_data.loc[all_ms1_data['file_category'].isin(['S1', 'ExCtrl'])]

    if logy:
        ms1_data[y_axis] = np.log10(ms1_data[y_axis])

    fig, axs = plt.subplots(2, sharex=True)
    title = compound_name

    pos_data = ms1_data.loc[ms1_data['polarity']=='POS']
    neg_data = ms1_data.loc[ms1_data['polarity']=='NEG']

    sample_pos_data = pos_data.loc[pos_data['file_category']=='S1']
    sample_neg_data = neg_data.loc[neg_data['file_category']=='S1']

    exctrl_pos_data = pos_data.loc[pos_data['file_category']=='ExCtrl']
    exctrl_neg_data = neg_data.loc[neg_data['file_category']=='ExCtrl']

    txctrl_pos_data = pos_data.loc[pos_data['file_category']=='TxCtrl']
    txctrl_neg_data = neg_data.loc[neg_data['file_category']=='TxCtrl']

    axs[0].scatter(sample_pos_data['run_num'], sample_pos_data[y_axis], color='orange', label='Sample')
    axs[0].scatter(exctrl_pos_data['run_num'], exctrl_pos_data[y_axis], color='gray', label='ExCtrl')
    axs[0].scatter(txctrl_pos_data['run_num'], txctrl_pos_data[y_axis], color='yellow', label='TxCtrl')
    
    axs[1].scatter(sample_neg_data['run_num'], sample_neg_data[y_axis], color='orange', label='Sample')
    axs[1].scatter(exctrl_neg_data['run_num'], exctrl_neg_data[y_axis], color='gray', label='ExCtrl')
    axs[1].scatter(txctrl_neg_data['run_num'], txctrl_neg_data[y_axis], color='yellow', label='TxCtrl')

    axs[0].margins(y=1)
    axs[1].margins(y=1)
    axs[0].set_title('POS')
    axs[1].set_title('NEG')

    axs[0].legend(loc="upper right", prop={'size': 6})
    
    fig.suptitle(title)
    fig.supxlabel('RUN_NUMBER')
    fig.supylabel(y_axis.upper())

    return fig

def _create_line_plots(ms1_tic_data, ax, cmap_name):

    n = ms1_tic_data.shape[0]
    cmap = cm.get_cmap(cmap_name, n+6)
    colors = cmap(np.arange(0,cmap.N)) 
    colors = colors[3:-3]

    ms1_tic_data.reset_index(drop=True, inplace=True)
    for idx, row in ms1_tic_data.iterrows():
        ax.plot(row.ms1_tic[1], row.ms1_tic[0], color=colors[idx], label='{group}_{run_num}'.format(group=row['group'], run_num=row['run_num']))

def _make_m1_tic_figure(all_ms1_tic_data, group_name, logy=False):
    ms1_tic_data = all_ms1_tic_data.loc[all_ms1_tic_data['file_category'].isin(['S1', 'ExCtrl'])]
    ms1_tic_data['ms1_tic'] = ms1_tic_data['ms1_tic'].apply(lambda x: np.asarray(eval(x)))

    if logy:
        ms1_tic_data['ms1_tic']= ms1_tic_data['ms1_tic'].apply(lambda x: np.array([np.log10(x[0]), x[1]]))

    fig, axs = plt.subplots(2, sharex=True)
    title = group_name

    pos_data = ms1_tic_data.loc[ms1_tic_data['polarity']=='POS']
    neg_data = ms1_tic_data.loc[ms1_tic_data['polarity']=='NEG']

    group_pos_data = pos_data.loc[pos_data['group']==group_name]
    group_neg_data = neg_data.loc[neg_data['group']==group_name]

    exctrl_pos_data = pos_data.loc[pos_data['file_category']=='ExCtrl']
    exctrl_neg_data = neg_data.loc[neg_data['file_category']=='ExCtrl']

    _create_line_plots(group_pos_data, axs[0], 'plasma')
    _create_line_plots(group_neg_data, axs[1], 'plasma')

    _create_line_plots(exctrl_pos_data, axs[0], 'Greys')
    _create_line_plots(exctrl_neg_data, axs[1], 'Greys')

    axs[0].margins(y=1)
    axs[1].margins(y=1)
    axs[0].set_title('POS')
    axs[1].set_title('NEG')

    axs[0].legend(loc="upper right", prop={'size': 6})
    
    fig.suptitle(title)
    fig.supxlabel('RETENTION_TIME')
    fig.supylabel('INTENSITY')

    return fig

def _label_outliers(data_df, y_axis, median_value, ax):
    """
    label points with run number outside of 50% or median value
    """
    for idx, row in data_df.iterrows():

        if y_axis == 'observed_intensity':
            if row[y_axis] > median_value * 10 or row[y_axis] < median_value * 0.1:
                ax.text(row.run_num, row[y_axis], row.run_num)
        if y_axis == 'ppm_error':
            if row[y_axis] > median_value + 1.75 or row[y_axis] < median_value - 2:
                ax.text(row.run_num, row[y_axis], row.run_num)
        if y_axis == 'retention_time':
            if row[y_axis] > median_value + 0.2 or row[y_axis] < median_value - 0.5:
                ax.text(row.run_num, row[y_axis], row.run_num)

        

#'observed_intensity' 'ppm_error' 'retention_time'

def _make_compound_ms1_figure(all_ms1_data, y_axis, compound_name, intensity_threshold=None, logy=False, label_outliers=True):

    ms1_data = all_ms1_data.loc[all_ms1_data['file_category'].isin(['ISTD', 'S1', 'ExCtrl'])]

    if logy:
        ms1_data[y_axis] = np.log10(ms1_data[y_axis])
        intensity_threshold = np.log10(intensity_threshold)

    fig, axs = plt.subplots(2, sharex=True)
    title = compound_name

    pos_data = ms1_data.loc[ms1_data['polarity']=='POS']
    neg_data = ms1_data.loc[ms1_data['polarity']=='NEG']

    istd_pos_data = pos_data.loc[pos_data['file_category']=='ISTD']
    istd_neg_data = neg_data.loc[neg_data['file_category']=='ISTD']

    sample_pos_data = pos_data.loc[pos_data['file_category']=='S1']
    sample_neg_data = neg_data.loc[neg_data['file_category']=='S1']

    exctrl_pos_data = pos_data.loc[pos_data['file_category']=='ExCtrl']
    exctrl_neg_data = neg_data.loc[neg_data['file_category']=='ExCtrl']

    txctrl_pos_data = pos_data.loc[pos_data['file_category']=='TxCtrl']
    txctrl_neg_data = neg_data.loc[neg_data['file_category']=='TxCtrl']

    sample_pos_y = np.array(sample_pos_data[y_axis].tolist())
    istd_pos_y = np.array(istd_pos_data[y_axis].tolist())
    exctrl_pos_y = np.array(exctrl_pos_data[y_axis].tolist())
    txctrl_pos_y = np.array(txctrl_pos_data[y_axis].tolist())

    sample_neg_y = np.array(sample_neg_data[y_axis].tolist())
    istd_neg_y = np.array(istd_neg_data[y_axis].tolist())
    exctrl_neg_y = np.array(exctrl_neg_data[y_axis].tolist())
    txctrl_neg_y = np.array(txctrl_neg_data[y_axis].tolist())

    if intensity_threshold is not None:
        sample_pos_col = np.where(sample_pos_y<intensity_threshold,'orangered', 'orange')
        istd_pos_col = np.where(istd_pos_y<intensity_threshold,'midnightblue', 'blue')
        exctrl_pos_col = np.where(exctrl_pos_y<intensity_threshold,'black', 'gray')
        txctrl_pos_col = np.where(txctrl_pos_y<intensity_threshold,'gold', 'yellow')

        sample_neg_col = np.where(sample_neg_y<intensity_threshold,'orangered', 'orange')
        istd_neg_col = np.where(istd_neg_y<intensity_threshold,'midnightblue', 'blue')
        exctrl_neg_col = np.where(exctrl_neg_y<intensity_threshold,'black', 'gray')
        txctrl_neg_col = np.where(txctrl_neg_y<intensity_threshold,'gold', 'yellow')
    else:
        sample_pos_col = 'orange'
        istd_pos_col = 'blue'
        exctrl_pos_col = 'gray'
        txctrl_pos_col = 'yellow'

        sample_neg_col = 'orange'
        istd_neg_col = 'blue'
        exctrl_neg_col = 'gray'
        txctrl_neg_col = 'yellow'

    axs[0].scatter(sample_pos_data['run_num'], sample_pos_y, color=sample_pos_col, label='Sample')
    axs[0].scatter(istd_pos_data['run_num'], istd_pos_y, color=istd_pos_col, label='ISTD')
    axs[0].scatter(exctrl_pos_data['run_num'], exctrl_pos_y, color=exctrl_pos_col, label='ExCtrl')
    axs[0].scatter(txctrl_pos_data['run_num'], txctrl_pos_y, color=txctrl_pos_col, label='TxCtrl')
    
    axs[1].scatter(sample_neg_data['run_num'], sample_neg_y, color=sample_neg_col, label='Sample')
    axs[1].scatter(istd_neg_data['run_num'], istd_neg_y, color=istd_neg_col, label='ISTD')
    axs[1].scatter(exctrl_neg_data['run_num'], exctrl_neg_y, color=exctrl_neg_col, label='ExCtrl')
    axs[1].scatter(txctrl_neg_data['run_num'], txctrl_neg_y, color=txctrl_neg_col, label='TxCtrl')

    if label_outliers is not None:
        pos_median_value = pos_data[y_axis].median()
        neg_median_value = neg_data[y_axis].median()
        _label_outliers(sample_pos_data, y_axis, pos_median_value, axs[0])
        _label_outliers(istd_pos_data, y_axis, pos_median_value, axs[0])
        _label_outliers(exctrl_pos_data, y_axis, pos_median_value, axs[0])
        _label_outliers(txctrl_pos_data, y_axis, pos_median_value, axs[0])

        _label_outliers(sample_neg_data, y_axis, neg_median_value, axs[1])
        _label_outliers(istd_neg_data, y_axis, neg_median_value, axs[1])
        _label_outliers(exctrl_neg_data, y_axis, neg_median_value, axs[1])
        _label_outliers(txctrl_neg_data, y_axis, neg_median_value, axs[1])

    axs[0].margins(y=1)
    axs[1].margins(y=1)
    axs[0].set_title('POS')
    axs[1].set_title('NEG')

    axs[0].hlines(y=np.mean(np.nan_to_num(istd_pos_y)), xmin=pos_data['run_num'].min(), xmax=pos_data['run_num'].max(), color='blue', linestyle=(0, (1, 10)), label='ISTD Mean')
    axs[1].hlines(y=np.mean(np.nan_to_num(istd_neg_y)), xmin=neg_data['run_num'].min(), xmax=neg_data['run_num'].max(), color='blue', linestyle=(0, (1, 10)), label='ISTD Mean')

    if intensity_threshold is not None:
        intensity_threshold = [intensity_threshold for i in range(len(ms1_data.run_num))]
        axs[0].plot(ms1_data['run_num'], intensity_threshold, color='red', linestyle='solid', label='Intensity Threshold')
        axs[1].plot(ms1_data['run_num'], intensity_threshold, color='red', linestyle='solid', label='Intensity Threshold')

    axs[0].legend(loc="upper right", prop={'size': 6})
    
    fig.suptitle(title)
    fig.supxlabel('RUN_NUMBER')
    fig.supylabel(y_axis.upper())

    return fig

def _make_compound_ms2_figure(ms2_data, y_axis, ms2_diagnostic, logy=False):
    
    if logy:
        ms2_data[y_axis] = np.log10(ms2_data[y_axis])

    pos_ms2_data = ms2_data.loc[ms2_data['polarity']=='POS']
    neg_ms2_data = ms2_data.loc[ms2_data['polarity']=='NEG']

    fig, axes = plt.subplots(2, sharex=True)
    title = 'Diagnostic MS2 Fragment Ions'

    for pos_frag_ion in ms2_diagnostic['POS']['diagnostic_ions']:
        pos_frag_ion_data = pos_ms2_data.loc[pos_ms2_data['theoretical_mz']==pos_frag_ion]
        axes[0].plot(pos_frag_ion_data['run_num'], pos_frag_ion_data[y_axis], label= pos_frag_ion, marker='o')

    for neg_frag_ion in ms2_diagnostic['NEG']['diagnostic_ions']:
        neg_frag_ion_data = neg_ms2_data.loc[neg_ms2_data['theoretical_mz']==neg_frag_ion]
        axes[1].plot(neg_frag_ion_data['run_num'], neg_frag_ion_data[y_axis], label= neg_frag_ion, marker='o')
        
    axes[0].margins(y=1)
    axes[1].margins(y=1)
    axes[0].set_title('POS')
    axes[1].set_title('NEG')
    axes[0].legend(loc="upper right", prop={'size': 6})
    axes[1].legend(loc="upper right", prop={'size': 6})

    fig.suptitle(title)
    fig.supxlabel('RUN_NUMBER')
    fig.supylabel(y_axis.upper())

    return fig


def make_ms1_qc_plots(ms1_data, atlas_df, qc_output_dir, logy=False):

    ms1_data.sort_values('run_num', inplace=True)

    if logy:
        logy_str = '_logy'
    else:
        logy_str = ''
    
    pdf = PdfPages(qc_output_dir + '/ms1_qc_plots{logy}.pdf'.format(logy=logy_str))

    for idx, row in atlas_df[atlas_df.signal_check==False].iterrows():

        compound_df = ms1_data.loc[ms1_data.compound_name == row.compound_name]

        fig1 = _make_compound_ms1_figure(compound_df, 'observed_intensity', row.compound_name, intensity_threshold=row.intensity_threshold, logy=logy)
        fig2 = _make_compound_ms1_figure(compound_df, 'ppm_error', row.compound_name)
        fig3 = _make_compound_ms1_figure(compound_df, 'retention_time', row.compound_name)

        pdf.savefig(fig1)
        pdf.savefig(fig2)
        pdf.savefig(fig3)

        plt.close('all')

    pdf.close()

def make_ms2_qc_plots(ms2_data, ms2_diagnostic, qc_output_dir, logy=False):

    ms2_data.sort_values('run_num', inplace=True)

    if logy:
        logy_str = '_logy'
    else:
        logy_str = ''
    
    pdf = PdfPages(qc_output_dir + '/ms2_qc_plots{logy}.pdf'.format(logy=logy_str))

    fig1 = _make_compound_ms2_figure(ms2_data, 'observed_intensity', ms2_diagnostic, logy=logy)
    fig2 = _make_compound_ms2_figure(ms2_data, 'ppm_error', ms2_diagnostic)

    pdf.savefig(fig1)
    pdf.savefig(fig2)

    plt.close('all')

    pdf.close()

def make_unlabeled_ms1_qc_plots(ms1_data, atlas_df, qc_output_dir, logy=False):

    ms1_data.sort_values('run_num', inplace=True)

    if logy:
        logy_str = '_logy'
    else:
        logy_str = ''
    
    pdf = PdfPages(qc_output_dir + '/ms1_unlabeled_intensity_plots{logy}.pdf'.format(logy=logy_str))

    for idx, row in atlas_df[atlas_df.signal_check==True].iterrows():

        compound_df = ms1_data.loc[ms1_data.compound_name == row.compound_name]

        fig1 = _make_unlabeled_m1_compound_figure(compound_df, 'observed_intensity', row.compound_name, logy=logy)
        fig2 = _make_unlabeled_m1_compound_figure(compound_df, 'ppm_error', row.compound_name)

        pdf.savefig(fig1)
        pdf.savefig(fig2)

        plt.close('all')

    pdf.close()

def make_ms1_tic_qc_plots(ms1_tic_data, qc_output_dir, logy=False, exclude_groups=['QC', 'ISTD', 'InjBl', 'ExCtrl']):

    if logy:
        logy_str = '_logy'
    else:
        logy_str = ''

    pdf = PdfPages(qc_output_dir + '/ms1_tic_plots{logy}.pdf'.format(logy=logy_str))

    groups = set(ms1_tic_data.group.tolist())
    [groups.discard(exclude_group) for exclude_group in exclude_groups]

    if len(groups) > 1:
        for group in groups:

            fig1 = _make_m1_tic_figure(ms1_tic_data, group, logy=logy)

            pdf.savefig(fig1)

            plt.close('all')
    
    else:
        print('Only one group detected, not generating overlay plots...')
        plt.close('all')
    
    pdf.close()