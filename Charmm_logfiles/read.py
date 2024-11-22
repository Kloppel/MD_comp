import os, re
import numpy as np 
import pandas as pd 

def read_charmm_blocks(logfile_path, entry_type):
    """
    Reads a CHARMM logfile and extracts blocks for the given entry type ('DYNA', 'AVER', or 'FLUC').

    Args:
        logfile_path (str): Path to the CHARMM logfile.
        entry_type (str): The type of entries to process ('DYNA', 'AVER', or 'FLUC').

    Returns:
        list: A list of block texts, each representing a block in the logfile.
    """
    first_phrase = f'{entry_type}>'
    blocks = []
    current_block = []
    collecting_block = False
    with open(logfile_path, 'r') as f:
        for line in f:
            if line.startswith(first_phrase):
                if collecting_block:
                    blocks.append(''.join(current_block))
                    current_block = [line]
                else:
                    collecting_block = True
                    current_block = [line]
            elif any(line.startswith(f'{entry_type} {suffix}>') for suffix in ['PROP', 'INTERN', 'CROSS', 'EXTERN', 'IMAGES', 'EWALD', 'CONSTR', 'PRESS']):
                if collecting_block:
                    current_block.append(line)
            else:
                if collecting_block:
                    current_block.append(line)
        if current_block:
            blocks.append(''.join(current_block))
    return blocks

def parse_block(block_text, entry_type):
    """
    Parses a block of text, mapping fixed column names to their corresponding values.

    Args:
        block_text (str): The block text to parse.
        entry_type (str): The type of entries ('DYNA', 'AVER', or 'FLUC').

    Returns:
        dict: A dictionary with variable names as keys and their corresponding values.
    """
    phrase_to_variables = {
        f'{entry_type}>': ['Step', 'Time', 'TOTEner', 'TOTKe', 'ENERgy', 'TEMPerature'],
        f'{entry_type} PROP>': ['GRMS', 'HFCTote', 'HFCKe', 'EHFCor', 'VIRKe'],
        f'{entry_type} INTERN>': ['BONDs', 'ANGLes', 'UREY-b', 'DIHEdrals', 'IMPRopers'],
        f'{entry_type} CROSS>': ['CMAPs', 'PMF1D', 'PMF2D', 'PRIMO'],
        f'{entry_type} EXTERN>': ['VDWaals', 'ELEC', 'HBONds', 'ASP', 'USER'],
        f'{entry_type} IMAGES>': ['IMNBvdw', 'IMELec', 'IMHBnd', 'RXNField', 'EXTElec'],
        f'{entry_type} EWALD>': ['EWKSum', 'EWSElf', 'EWEXcl', 'EWQCor', 'EWUTil'],
        f'{entry_type} CONSTR>': ['HARMonic', 'CDIHedral', 'CIC', 'RESDistance', 'NOE'],
        f'{entry_type} PRESS>': ['VIRE', 'VIRI', 'PRESSE', 'PRESSI', 'VOLUme']
    }
    data = {}
    lines = block_text.strip().split('\n')
    for line in lines:
        for phrase, variables in phrase_to_variables.items():
            if line.startswith(phrase):
                line_content = line[len(phrase):].strip()
                values = re.findall(r'[-+]?\d*\.\d+|\d+', line_content)
                if len(values) != len(variables):
                    print(f"Warning: Number of values ({len(values)}) does not match number of variables ({len(variables)}) for line: {line}")
                    min_len = min(len(values), len(variables))
                    values = values[:min_len]
                    variables = variables[:min_len]
                for var, val in zip(variables, values):
                    try:
                        data[var] = float(val)
                    except ValueError:
                        data[var] = val
                break
    return data

def process_charmm_logfile(logfile_path, entry_type):
    """
    Processes a CHARMM logfile, extracts blocks of a given type (DYNA, AVER, or FLUC),
    parses them using parse_block, and returns a DataFrame.

    Args:
        logfile_path (str): Path to the CHARMM logfile.
        entry_type (str): The type of entries to process ('DYNA', 'AVER', or 'FLUC').

    Returns:
        DataFrame: A pandas DataFrame containing the parsed data.
    """
    blocks = read_charmm_blocks(logfile_path, entry_type)
    data_list = []
    for block_text in blocks:
        data = parse_block(block_text, entry_type)
        if data:
            data_list.append(data)
    df = pd.DataFrame(data_list)
    if 'Step' in df.columns:
        df = df.sort_values(by='Step').reset_index(drop=True)
    df = df.drop_duplicates()
    return df

def charmm_logfile_all(logfile_path):
    #logfile_path = '1eru_eq.logfile'
    df_dyna = process_charmm_logfile(logfile_path, 'DYNA')
    df_aver = process_charmm_logfile(logfile_path, 'AVER')
    df_fluc = process_charmm_logfile(logfile_path, 'FLUC')
    return df_dyna, df_aver, df_fluc