import numpy as np
import fileinput
import re
import os
import pandas as pd
import inspect
import types
from typing import cast

def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def validate_shell_in_orca_main_basis(s):
    try:
        return int(s)
    except (ValueError,TypeError):
        return s

# Read the provided orbital basis set file as input and strip any whitespaces
def read_and_format_orca_basis_set(original_basis):
    formatted_basis = '%s-formatted_basis' % original_basis
    with open('%s' % original_basis, "r") as infile, open('%s-formatted_basis' % original_basis, "w") as outfile:
        for line in fileinput.input(original_basis):
            line = re.sub(' +', ',', line)
            outfile.write(''.join(line))
    # Remove $DATA and $END labels appearing in the ORCA basis set file
    with open('%s-formatted_basis' % original_basis, "r") as infile2:
        lines = infile2.readlines()
    with open('%s-formatted_basis' % original_basis, "w") as infile2:
        for line in lines:
            line = line.replace('$DATA','')
            line = line.replace('$END','')
            if re.match(r'^\s*$', line):
                pass
            else:
                infile2.write(line)
    # Remove now the first row that contains the atom name, e.g, ALUMINIUM
    with open('%s-formatted_basis' % original_basis,'r') as fin:
        data = fin.read().splitlines(True)
    with open('%s-formatted_basis' % original_basis,'w') as out:
        out.writelines(data[1:])

    # Ensure that all shells are assigned to letters in the first column
    # of the main basis set, thus replace any 'int' by the letter of the corresponding shell.
    # Also, save the column with shell definitions in a single csv file
    current_directory = (os.getcwd())
    main_basis = pd.read_csv(f'{current_directory}/{formatted_basis}', sep=',',skiprows=0,names=['Shell','Orbital_Exponent','Contract_Coeff'])
    shell_basis = main_basis['Shell']
    shell_basis_formatted = [validate_shell_in_orca_main_basis(s) for s in shell_basis]
    for i, n in enumerate(shell_basis_formatted):
        if type(n) == int:
           shell_basis_formatted[i] = shell_basis_formatted[i-1]
    shell_basis_formatted = pd.DataFrame(shell_basis_formatted)
    # Remove older single column with shells and replace with the generated one
    main_basis_formatted = pd.concat([shell_basis_formatted, main_basis['Orbital_Exponent'], main_basis['Contract_Coeff']],axis=1)
    # Drop NaN values
    main_basis_formatted = main_basis_formatted.dropna()
    # Save the formatted basis
    main_basis_formatted.to_csv(formatted_basis,index=False,header=['Shell','Orbital_Exponent','Contraction_Coeff'])

#################################################################################
# Read the provided orbital basis set gbs file as input and strip any whitespaces
def read_and_format_gbs_basis_set(original_basis):
    formatted_basis = '%s-formatted_basis' % original_basis
    with open('%s' % original_basis, "r") as infile, open('%s-formatted_basis' % original_basis, "w") as outfile:
        for line in fileinput.input(original_basis):
            line = re.sub(' +', ',', line)
            line = re.sub("S.*","S",line)
            line = re.sub("P.*","P",line)
            line = re.sub("D.*","D",line)
            line = re.sub("F.*","F",line)
            line = re.sub("G.*","G",line)
            line = re.sub("H.*","H",line)
            outfile.write(''.join(line))
    # Remove asterisks in the gbs basis set file
    with open('%s-formatted_basis' % original_basis, "r") as infile2:
        lines = infile2.readlines()
    with open('%s-formatted_basis' % original_basis, "w") as infile2:
        for line in lines:
            line = line.replace('****','')
            if re.match(r'^\s*$', line):
                pass
            else:
                infile2.write(line)
    # Remove now the first row that contains the atom name, e.g, ALUMINIUM
    with open('%s-formatted_basis' % original_basis,'r') as fin:
        data = fin.read().splitlines(True)
    with open('%s-formatted_basis' % original_basis,'w') as out:
        out.writelines(data[1:])
    # Ensure that all shells are assigned to letters in the first column
    # of the main basis set, thus replace any 'int' by the letter of the corresponding shell.
    # Also, save the column with shell definitions in a single csv file
    current_directory = (os.getcwd())
    main_basis = pd.read_csv(f'{current_directory}/{formatted_basis}', sep=',',skiprows=0,names=['Shell','Orbital_Exponent','Contract_Coeff'])
    shell_basis = main_basis['Shell']
    shell_basis = shell_basis.fillna(1)
    shell_basis_formatted = [validate_shell_in_orca_main_basis(s) for s in shell_basis]
    for i, n in enumerate(shell_basis_formatted):
        if type(n) == int:
           shell_basis_formatted[i] = shell_basis_formatted[i-1]
    shell_basis_formatted = pd.DataFrame(shell_basis_formatted)
    # Remove older single column with shells and replace with the generated one
    main_basis_formatted = pd.concat([shell_basis_formatted, main_basis['Orbital_Exponent'], main_basis['Contract_Coeff']],axis=1)
    # Drop NaN value
    main_basis_formatted = main_basis_formatted.dropna()
    # Save the formatted basis
    main_basis_formatted.to_csv(formatted_basis,index=False,header=['Shell','Orbital_Exponent','Contraction_Coeff'])


def concatenate_and_short(a, b, kind='mergesort'):
    c = np.concatenate((a, b))
    c.sort(kind=kind)
    flag = np.ones(len(c),dtype=bool)
    np.not_equal(c[1:],c[:-1],out=flag[1:])
    return c[flag]

def write_exponents_summary_to_gbs_file(exponents_summary,exponents_summary_name,SHELL_NAME):
    # Save all exponents by shell type in descending order:
    if exponents_summary.size != 0:
        with open(f'{exponents_summary_name}','w') as f:
            f.write('%s\n' % SHELL_NAME)
            for item in (np.sort(exponents_summary)[::-1]):
                f.write(item.astype(str))
                f.write('\n')
    # Combine the generated AutoCABS exponents into a single intermediate gbs file
    if os.path.exists(f'{exponents_summary_name}') == True:
        with open(f'{exponents_summary_name}','r') as f, open(f'{exponents_summary_name}_gbs','w') as out:
            header = f.readline()
            for i, line in enumerate(f):
                line = line.rstrip('\n')
                header = header.rstrip('\n')
                out.write('%s  1  1.00\n' % header)
                out.write('    %s \t 1.0000000\n' % line)
        os.remove(exponents_summary_name)

def write_exponents_summary_to_orca_file(exponents_summary,exponents_summary_name,SHELL_NAME):
    # Save all exponents by shell type in descending order:
    if exponents_summary.size != 0:
        with open(f'{exponents_summary_name}','w') as f:
            f.write('%s\n' % SHELL_NAME)
            for item in (np.sort(exponents_summary)[::-1]):
                f.write(item.astype(str))
                f.write('\n')
    # Combine the generated AutoCABS exponents into a single Orca-compatible basis set file
    if os.path.exists(f'{exponents_summary_name}') == True:
        with open(f'{exponents_summary_name}','r') as f, open(f'{exponents_summary_name}_gbs','w') as out:
            header = f.readline()
            for i, line in enumerate(f):
                line = line.rstrip('\n')
                header = header.rstrip('\n')
                out.write('%s  1\n' % header)
                out.write('1    %s  1.0000000\t \n' % line)
        os.remove(exponents_summary_name)


def write_exponents_summary_to_tm_file(exponents_summary,exponents_summary_name,SHELL_NAME):
    # Save all exponents by shell type in descending order:
    if exponents_summary.size != 0:
        with open(f'{exponents_summary_name}','w') as f:
            f.write('%s\n' % SHELL_NAME)
            for item in (np.sort(exponents_summary)[::-1]):
                f.write(item.astype(str))
                f.write('\n')
    # Combine the generated AutoCABS exponents into a single Turbomole-compatible basis set file
    if os.path.exists(f'{exponents_summary_name}') == True:
        with open(f'{exponents_summary_name}','r') as f, open(f'{exponents_summary_name}_tm','w') as out:
            header = f.readline()
            for i, line in enumerate(f):
                line = line.rstrip('\n')
                header = header.rstrip('\n')
                out.write('   1  %s\n' % header)
                out.write('     %s  1.0000000\t \n' % line)
        os.remove(exponents_summary_name)
