#!/usr/bin/env python
import json
import os
import re
import warnings

import matplotlib.pyplot as plt
import neurokit as nk
import numpy as np
import pandas as pd


# functions to analyze gsr data
def analyze_gsr(acqs, txts, outdir):
    """
    Main function to analyze the gsr data
    """
    group_tsv = os.path.join(outdir, "group.tsv")
    acq_nonmatches = []
    group_msr = []
    for acq in acqs:
        acq_txt_dict, acq_nonmatch, txts = match_acq_txt(acq, txts)
        if acq_txt_dict is None:
            acq_nonmatches.append(acq_nonmatch)
            continue
        
        fname_dict = gen_filenames(acq_txt_dict)
        
        if (
            os.path.isfile(fname_dict['fig']) and
            os.path.isfile(fname_dict['data']) and
            os.path.isfile(fname_dict['smry'])
           ):
            with open(fname_dict['smry']) as json_file:  
                sub_dict = json.load(json_file)
        else:
            acq_txt_dict_proc = filter_gsr(acq_txt_dict)
            sub_dict = write_results(acq_txt_dict_proc, fname_dict)
        
        group_msr.append(sub_dict)
  
    group_df = pd.DataFrame(group_msr)
    group_df.to_csv(group_tsv, sep="\t")
    txt_nonmatches = txts
    
    return acq_txt_dict_proc, acq_nonmatches, txt_nonmatches, group_df


def match_acq_txt(acq, txts):
    """
    Match the unclearly named acqknowledge file with
    the clearly named txt file
    """
    # minimum length of the data file (data for 1 second)
    min_len = 200
    # connects the acqknowledge file to the txt file
    acq_txt_dict = {}
    # read the acqknowledge file
    acq_df, samp_rate = nk.bio.bio_data.read_acqknowledge(acq)
    # there should be enough rows in the dataframes
    if acq_df.shape[0] < min_len:
        return None, acq, txts
            
    # big loop to see which txt file goes with the acknowledge file
    for txt in txts:
        # only look at txt files in the same directory as acknowledge file
        if os.path.dirname(acq) != os.path.dirname(txt):
            continue
        # attempt 1 to read the txt file
        # assumes no header information and no time column
        txt_df = pd.read_csv(txt,
        header=None,
        index_col=False,
        names=acq_df.columns,
        sep="\t")

        # the file may have a header
        if np.any(txt_df.iloc[0,:].isnull()):
            # assume the txt file has a header
            try:
                txt_df = pd.read_csv(txt,
                    header=14,
                    index_col=0,
                    delim_whitespace=True,
                    names=acq_df.columns)
            except:
                # assume the acknowledge file does not have the "Scanner TTl column"
                columns = list(acq_df.columns)
                columns.insert(-1, "Scanner TTL")
                txt_df = pd.read_csv(txt,
                    header=14,
                    index_col=0,
                    delim_whitespace=True,
                    names=columns)
            # assume in the two above cases,
            # there is a time column being treated as the index
            txt_df.reset_index(drop=True, inplace=True)
        # select this column for comparisons
        # assume the txt and acq have at least min_len rows
        a = acq_df['PPG100C'].reset_index(drop=True)[:min_len]
        t = txt_df['PPG100C'][:min_len]

        # the columns cannot be identical due to rounding error,
        # but can be close
        if np.all(np.isclose(a, t)):
            acq_txt_dict[acq] = {'txt': txt, 'df': acq_df}
            txts.remove(txt)
            return acq_txt_dict, None, txts
    
    # nothing matches
    return None, acq, txts


def filter_gsr(acq_txt_dict):
    """
    Filter/Preprocess the GSR data
    """
    key = list(acq_txt_dict.keys())[0]
    
    acq_df = acq_txt_dict[key]['df']
    # where the Scanner Trigger is 0 initiates a slice acquition in the scanner
    # this may be false in some scenerios (e.g. when scan trigger is always 0)
    start = acq_df["Scanner Trigger"].where(acq_df["Scanner Trigger"]==0.0).first_valid_index()
    end = acq_df["Scanner Trigger"].where(acq_df["Scanner Trigger"]==0.0).last_valid_index()
    if start is None or end is None:
        warnings.warn("Scan Trigger Not Recognized, keeping all rows")
        acq_cut_df = acq_df
    else:
        acq_cut_df = acq_df.iloc[(acq_df.index >= start) & (acq_df.index <= end)]

    # process the GSR data
    res = nk.bio_eda.eda_process(acq_df['GSR100C'], filter_type='butter',
                                 band="lowpass", order=1, frequency=1,
                                 sampling_rate=200)

    acq_txt_dict[key]['df'] = res['df']
    return acq_txt_dict


def gen_filenames(acq_txt_dict):
    """
    Generate the output file names
    """
    # pattern for matching atrain files
    atrain_ptrn = re.compile(r".*/"
                             r"(?P<sub_id>sub-[A-Za-z0-9]+)"
                             r"/(?P<ses_id>ses-[A-Za-z0-9]+)"
                             r"/(?P<task_id>[A-Za-z]+)[1-9]"
                             r"(_(?P<run_id>[0-9]))?.txt")
    # pattern for matching extend files
    extend_ptrn = re.compile(r".*/"
                             r"(?P<sub_id>sub-[A-Za-z0-9]+)"
                             r"(_(?P<ses_id>ses-[A-Za-z0-9]+))"
                             r"(_(?P<task_id>task-[A-Za-z0-9]+))"
                             r"(_(?P<run_id>run-[0-9]+))?"
                             r"(_physio)?.txt")
    for acq, tgt_dct in acq_txt_dict.items():
        atrain_mch = atrain_ptrn.match(tgt_dct['txt'])
        extend_mch = extend_ptrn.match(tgt_dct['txt'])
        
        if atrain_mch:
            fdict = atrain_mch.groupdict()
            # make task lowercase
            fdict['task_id'] = fdict['task_id'].lower()
            # make the runs 1 or 2 (since there are only 2 runs per session)
            if fdict['run_id'] == "3" or fdict['run_id'] == "4":
                fdict['run_id'] = int(fdict['run_id']) - 2
        elif extend_mch:
            fdict = extend_mch.groupdict()
            # make run an integer if necessary
            if fdict["run_id"]:
                fdict["run_id"] = int(fdict["run_id"].lstrip("run-"))
        else:
            raise("FileName did not match either atrain or extend pattern")

        # strip the keys from the labels
        fdict["sub_id"] = fdict["sub_id"].lstrip("sub-")
        fdict["ses_id"] = fdict["ses_id"].lstrip("ses-")
        fdict["task_id"] = fdict["task_id"].lstrip("task-")
        
        # template output file changes depending if run is a key or not
        if not fdict["run_id"]:
            tmplt = os.path.join(outdir,
                                 "sub-{sub_id}",
                                 "ses-{ses_id}",
                                 "sub-{sub_id}_ses-{ses_id}_task-{task_id}_{typ}.{ext}")
        elif fdict["run_id"]:
            tmplt = os.path.join(outdir,
                                 "sub-{sub_id}",
                                 "ses-{ses_id}",
                                 "sub-{sub_id}_ses-{ses_id}_task-{task_id}_run-{run_id}_{typ}.{ext}")
        fig_file = tmplt.format(**fdict, typ="qa", ext="svg")
        data_file = tmplt.format(**fdict, typ="physio", ext="tsv")
        sum_file = tmplt.format(**fdict, typ="summary", ext="json")
        
        return {'fig': fig_file, 'data': data_file, 'smry': sum_file, 'fdict': fdict}


def write_results(acq_txt_dict, filenames):
    """
    write out the data from the processed GSR data
    """
    # make directory to place results
    os.makedirs(os.path.dirname(filenames['data']), exist_ok=True)
    key = list(acq_txt_dict.keys())[0]
    
    fdict = filenames['fdict']
    out_df = acq_txt_dict[key]['df']
    out_df.to_csv(filenames['data'])
    title = os.path.basename(filenames['fig'].split('.')[0])
    out_df.plot(y=["EDA_Raw", "EDA_Filtered"], title=title)
    plt.savefig(filenames['fig'])
    plt.clf()
    plt.close()

    sub_dict = {
        "participant_id": fdict["sub_id"],
        "session_id": fdict["ses_id"],
        "task_id": fdict["task_id"],
        "run_id": fdict.get("run_id"),
        "EDA_mean": out_df["EDA_Filtered"].mean(),
        "EDA_median": out_df["EDA_Filtered"].median(),
        "EDA_std": out_df["EDA_Filtered"].std()
    }

    with open(filenames['smry'], 'w') as outfile:  
        json.dump(sub_dict, outfile)
    
    return sub_dict

        
if __name__ == '__main__':
    from argparse import ArgumentParser
    from glob import glob

    parser = ArgumentParser(description='physio.py: analyze physio data')

    parser.add_argument('indir', action='store',
                        help='path to the top level of the data directory')
    parser.add_argument('outdir', action='store',
                        help='path to output directory')
    
    opts = parser.parse_args()

    acqs = glob(os.path.join(opts.indir, "sub-*", "**", "*.acq"))
    txts = glob(os.path.join(opts.indir, "sub-*", "**", "*.txt"))
    outdir = os.path.abspath(opts.outdir)

    analyze_gsr(acqs, txts, outdir)
