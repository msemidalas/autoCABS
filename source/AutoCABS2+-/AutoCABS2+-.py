import sys
import csv
import os
import glob
from collections import defaultdict
import numpy as np
import pandas as pd
import basis_set_exchange as bse
from myfunctions import concatenate_and_short
from myfunctions import read_and_format_orca_basis_set
from myfunctions import read_and_format_gbs_basis_set
from myfunctions import unique_rows
from myfunctions import write_exponents_summary_to_gbs_file
from myfunctions import write_exponents_summary_to_orca_file
from myfunctions import write_exponents_summary_to_tm_file
from chemtools.basisset import BasisSet, merge
# Requires MolSSI-BSE (Basis set exchange) version 0.9 and Chemtools version 0.9.2

# Read and edit one MOLPRO (e.g., c.mpro) or ORCA basis (e.g.,c.orca).
# Remove any generated AutoCABS files that might be here from a previous run.
for f in glob.glob("*-AutoCABS*"):
    os.remove(f)

# Locate the MOLPRO or ORCA basis set files
molpro_filenames = glob.glob('**/*.mpro',recursive=True)
filename_molpro_original = len(molpro_filenames)
orca_filenames = glob.glob('**/*.orca',recursive=True)
filename_orca_original = len(orca_filenames)
if filename_molpro_original!=0:
    print('Basis set file in a MOLPRO format (.mpro) was found.')
    element_ID = (molpro_filenames[0]).split(".")[0].title()
    print('Chosen element is: %s' % element_ID)
    formatted_basis_mpro = BasisSet.from_file('%s' % molpro_filenames[0], fmt='molpro', name='user-def-basis')
    formatted_basis_gbs = formatted_basis_mpro.to_gaussian()
    with open('%s.gbs' % element_ID,"w") as outfile:
        outfile.write(formatted_basis_gbs)
    read_and_format_gbs_basis_set('%s.gbs' % element_ID)
    formatted_basis = '%s.gbs-formatted_basis' % element_ID
elif filename_orca_original!=0:
    print('Basis set file in an ORCA format was found.')
    element_ID = (orca_filenames[0]).split(".")[0].title()
    print('Chosen element is:%s' % element_ID)
    formatted_basis = '%s-formatted_basis' % orca_filenames[0]
    read_and_format_orca_basis_set(orca_filenames[0])
elif (filename_molpro_original==0 and filename_orca_original==0):
    print('No basis set file was found.')
    sys.exit()

# Store all shells, exponents, and contraction coefficients in a dictionary of lists
mydictionary = defaultdict(list)
with open('%s' % formatted_basis) as in_file:
    csv_reader = csv.DictReader(in_file, delimiter=",")
    for row in csv_reader:
        for key, values in row.items():
            mydictionary[key].append(values)

df = pd.DataFrame({'Shell': mydictionary['Shell'], 'Orbital_Exponent': pd.to_numeric(mydictionary['Orbital_Exponent']),
                   'Contraction_Coeff': pd.to_numeric(mydictionary['Contraction_Coeff'])})
df['S_Type_Exp'] = np.where(df['Shell'] == 'S',1,0)
df['P_Type_Exp'] = np.where(df['Shell'] == 'P',1,0)
df['D_Type_Exp'] = np.where(df['Shell'] == 'D',1,0)
df['F_Type_Exp'] = np.where(df['Shell'] == 'F',1,0)
df['G_Type_Exp'] = np.where(df['Shell'] == 'G',1,0)
df['H_Type_Exp'] = np.where(df['Shell'] == 'H',1,0)
df['I_Type_Exp'] = np.where(df['Shell'] == 'I',1,0)
df['Contraction'] = np.where(df['Contraction_Coeff'] == 1,1,0)
pd.set_option('display.max_columns',None)

# Group exponents in dataframes by shells.
S_orbital_exponents = df[df['S_Type_Exp']==1]
P_orbital_exponents = df[df['P_Type_Exp']==1]
D_orbital_exponents = df[df['D_Type_Exp']==1]
F_orbital_exponents = df[df['F_Type_Exp']==1]
G_orbital_exponents = df[df['G_Type_Exp']==1]
H_orbital_exponents = df[df['H_Type_Exp']==1]
I_orbital_exponents = df[df['I_Type_Exp']==1]

# Group exponents based on whether their contraction coeffs. are equal to one or not.
def collect_exponents_by_contraction(orbital_exponents_df):
    with_contraction_coeff_1 = orbital_exponents_df[orbital_exponents_df['Contraction'] == 1]
    with_contraction_coeff_1_np = np.array(with_contraction_coeff_1['Orbital_Exponent'], dtype="float64")
    if with_contraction_coeff_1_np.size != 0:
        with_contraction_coeff_1_np = np.sort(unique_rows(np.vstack(with_contraction_coeff_1_np)), axis=-1)

    with_contraction_coeff_not_1 = orbital_exponents_df[orbital_exponents_df['Contraction'] != 1]
    with_contraction_coeff_not_1_np = np.array(with_contraction_coeff_not_1['Orbital_Exponent'], dtype="float64")
    if with_contraction_coeff_not_1_np.size != 0:
        with_contraction_coeff_not_1_np = np.sort(unique_rows(np.vstack(with_contraction_coeff_not_1_np)), axis=-1)

    return with_contraction_coeff_1_np, with_contraction_coeff_not_1_np

# S exponents
S_orbital_exponents_with_contraction_coeff_1_np, S_orbital_exponents_with_contraction_coeff_not_1_np = collect_exponents_by_contraction(S_orbital_exponents)
print('\n S main orbital exponents with contraction coefficient equal to 1:\n', S_orbital_exponents_with_contraction_coeff_1_np)
print('\n S main orbital exponents with contraction coefficient not equal to 1:\n', S_orbital_exponents_with_contraction_coeff_not_1_np)

# P exponents
P_orbital_exponents_with_contraction_coeff_1_np, P_orbital_exponents_with_contraction_coeff_not_1_np = collect_exponents_by_contraction(P_orbital_exponents)
print('\n P main orbital exponents with contraction coefficient equal to 1:\n', P_orbital_exponents_with_contraction_coeff_1_np)
print('\n P main orbital exponents with contraction coefficient not equal to 1:\n', P_orbital_exponents_with_contraction_coeff_not_1_np)

# D exponents
D_orbital_exponents_with_contraction_coeff_1_np, D_orbital_exponents_with_contraction_coeff_not_1_np = collect_exponents_by_contraction(D_orbital_exponents)
print('\n D main orbital exponents with contraction coefficient equal to 1:\n', D_orbital_exponents_with_contraction_coeff_1_np)
print('\n D main orbital exponents with contraction coefficient not equal to 1:\n', D_orbital_exponents_with_contraction_coeff_not_1_np)

# F exponents
F_orbital_exponents_with_contraction_coeff_1_np, F_orbital_exponents_with_contraction_coeff_not_1_np = collect_exponents_by_contraction(F_orbital_exponents)
print('\n F main orbital exponents with contraction coefficient equal to 1:\n', F_orbital_exponents_with_contraction_coeff_1_np)
print('\n F main orbital exponents with contraction coefficient not equal to 1:\n', F_orbital_exponents_with_contraction_coeff_not_1_np)

# G exponents
G_orbital_exponents_with_contraction_coeff_1_np, G_orbital_exponents_with_contraction_coeff_not_1_np = collect_exponents_by_contraction(G_orbital_exponents)
print('\n G main orbital exponents with contraction coefficient equal to 1:\n', G_orbital_exponents_with_contraction_coeff_1_np)
print('\n G main orbital exponents with contraction coefficient not equal to 1:\n', G_orbital_exponents_with_contraction_coeff_not_1_np)

# H exponents
H_orbital_exponents_with_contraction_coeff_1_np, H_orbital_exponents_with_contraction_coeff_not_1_np = collect_exponents_by_contraction(H_orbital_exponents)
print('\n H main orbital exponents with contraction coefficient equal to 1:\n', H_orbital_exponents_with_contraction_coeff_1_np)
print('\n H main orbital exponents with contraction coefficient not equal to 1:\n', H_orbital_exponents_with_contraction_coeff_not_1_np)

# I exponents
I_orbital_exponents_with_contraction_coeff_1_np, I_orbital_exponents_with_contraction_coeff_not_1_np = collect_exponents_by_contraction(I_orbital_exponents)                                            
print('\n I main orbital exponents with contraction coefficient equal to 1:\n', I_orbital_exponents_with_contraction_coeff_1_np)                                                                        
print('\n I main orbital exponents with contraction coefficient not equal to 1:\n', I_orbital_exponents_with_contraction_coeff_not_1_np)

#############################################################################################################################
# Check for general contraction; extract exponents from all primitives per shell.
# For segmented contraction, obtain exponents with a coefficient of 1 and the outermost exponent with a non-1 coefficient.
#############################################################################################################################
def process_primitives(with_contraction_coeff_1, with_contraction_coeff_not_1):
    if with_contraction_coeff_1.size == 0:
        return with_contraction_coeff_not_1
    else:
        outermost_decontracted = np.setdiff1d(with_contraction_coeff_not_1.flatten(), with_contraction_coeff_1.flatten())
        if outermost_decontracted.size != 0:
            return concatenate_and_short(with_contraction_coeff_1.flat, outermost_decontracted[0,].flat, kind='mergesort')
        else:
            return with_contraction_coeff_1

# Collect S, P, D, F, G, H, and I exponents
S_orbital_exponents_np = process_primitives(S_orbital_exponents_with_contraction_coeff_1_np, S_orbital_exponents_with_contraction_coeff_not_1_np)
P_orbital_exponents_np = process_primitives(P_orbital_exponents_with_contraction_coeff_1_np, P_orbital_exponents_with_contraction_coeff_not_1_np)
D_orbital_exponents_np = process_primitives(D_orbital_exponents_with_contraction_coeff_1_np, D_orbital_exponents_with_contraction_coeff_not_1_np)
F_orbital_exponents_np = process_primitives(F_orbital_exponents_with_contraction_coeff_1_np, F_orbital_exponents_with_contraction_coeff_not_1_np)
G_orbital_exponents_np = process_primitives(G_orbital_exponents_with_contraction_coeff_1_np, G_orbital_exponents_with_contraction_coeff_not_1_np)
H_orbital_exponents_np = process_primitives(H_orbital_exponents_with_contraction_coeff_1_np, H_orbital_exponents_with_contraction_coeff_not_1_np)
I_orbital_exponents_np = process_primitives(I_orbital_exponents_with_contraction_coeff_1_np, I_orbital_exponents_with_contraction_coeff_not_1_np)
##############################################################################
# If outer shell contains a single main orbital exponent, then discard it.
# Instead, multiply by 1.5 the L-1 orbital exponents, and assign them to the outer shell:
if D_orbital_exponents_np.size == 1:
  D_orbital_exponents_np = 1.5 * P_orbital_exponents_np
  print('\nD orbital exponents were taken as multiples of P-type shells.\n')
if F_orbital_exponents_np.size == 1:
  F_orbital_exponents_np = 1.5 * D_orbital_exponents_np
  print('\nF orbital exponents were taken as multiples of D-type shells.\n')
if G_orbital_exponents_np.size == 1:
  G_orbital_exponents_np = 1.5 * F_orbital_exponents_np
  print('\nG orbital exponents were taken as multiples of F-type shells.\n')
if H_orbital_exponents_np.size == 1:
  H_orbital_exponents_np = 1.5 * G_orbital_exponents_np
  print('\nH orbital exponents were taken as multiples of G-type shells.\n')
if I_orbital_exponents_np.size == 1:                                                                                                                                                                    
  I_orbital_exponents_np = 1.5 * H_orbital_exponents_np                                                                                                                                                 
  print('\nI orbital exponents were taken as multiples of H-type shells.\n')
#################################################################################

######################################################
# Geometric mean calculation of main orbital exponents
######################################################
if S_orbital_exponents_np.size not in {0,1}:
    geom_mean_S = (S_orbital_exponents_np[:-1] * S_orbital_exponents_np[1:]) ** (1 / 2)
    geom_mean_S = unique_rows(np.vstack(geom_mean_S))
    print('\nGeometric mean of S orbital exponents\n', geom_mean_S)

if P_orbital_exponents_np.size not in {0,1}:
    geom_mean_P = (P_orbital_exponents_np[:-1] * P_orbital_exponents_np[1:]) ** (1 / 2)
    geom_mean_P = unique_rows(np.vstack(geom_mean_P))
    print('\nGeometric mean of P orbital exponents\n', geom_mean_P)
else:
    geom_mean_P = np.empty(0)
    print('\nGeometric mean of P orbital exponents cannot be calculated from current main orbital exponents:', geom_mean_P)

if D_orbital_exponents_np.size not in {0,1}:
    geom_mean_D = (D_orbital_exponents_np[:-1] * D_orbital_exponents_np[1:]) ** (1 / 2)
    geom_mean_D = unique_rows(np.vstack(geom_mean_D))
    print('\nGeometric mean of D orbital exponents\n', geom_mean_D)
else:
    geom_mean_D = np.empty(0)
    print('\nGeometric mean of D orbital exponents cannot be calculated from current main orbital exponents:', geom_mean_D)

if F_orbital_exponents_np.size not in {0,1}:
    geom_mean_F = (F_orbital_exponents_np[:-1] * F_orbital_exponents_np[1:]) ** (1 / 2)
    geom_mean_F = unique_rows(np.vstack(geom_mean_F))
    print('\nGeometric mean of F orbital exponents\n', geom_mean_F)
else:
    geom_mean_F = np.empty(0)
print('\nGeometric mean of F orbital exponents cannot be calculated from current main orbital exponents:', geom_mean_F)

if G_orbital_exponents_np.size not in {0,1}:
    geom_mean_G = (G_orbital_exponents_np[:-1] * G_orbital_exponents_np[1:]) ** (1 / 2)
    geom_mean_G = unique_rows(np.vstack(geom_mean_G))
    print('\nGeometric mean of G orbital exponents\n', geom_mean_G)
else:
    geom_mean_G = np.empty(0)
    print('\nGeometric mean of G orbital exponents cannot be calculated from current main orbital exponents:', geom_mean_G)

if H_orbital_exponents_np.size not in {0,1}:
    geom_mean_H = (H_orbital_exponents_np[:-1] * H_orbital_exponents_np[1:]) ** (1 / 2)
    geom_mean_H = unique_rows(np.vstack(geom_mean_H))
    print('\nGeometric mean of H orbital exponents\n', geom_mean_H)
else:
    geom_mean_H = np.empty(0)
    print('\nGeometric mean of H orbital exponents cannot be calculated from current main orbital exponents:', geom_mean_H)

if I_orbital_exponents_np.size not in {0,1}:                                                                                                                                                            
    geom_mean_I = (I_orbital_exponents_np[:-1] * I_orbital_exponents_np[1:]) ** (1 / 2)                                                                                                                 
    geom_mean_I = unique_rows(np.vstack(geom_mean_I))                                                                                                                                                   
    print('\nGeometric mean of I orbital exponents\n', geom_mean_I)                                                                                                                                     
else:                                                                                                                                                                                                   
    geom_mean_I = np.empty(0)                                                                                                                                                                           
    print('\nGeometric mean of I orbital exponents cannot be calculated from current main orbital exponents:', geom_mean_I)

# Summary with generated exponents
print('\nSummary of generated AutoCABS exponents')
print('S\n', np.sort(geom_mean_S)[::-1])
print('P\n', np.sort(geom_mean_P)[::-1])
print('D\n', np.sort(geom_mean_D)[::-1])
print('F\n', np.sort(geom_mean_F)[::-1])
print('G\n', np.sort(geom_mean_G)[::-1])
print('H\n', np.sort(geom_mean_H)[::-1])
print('I\n', np.sort(geom_mean_I)[::-1])
# ##############################################################################################
# # Addition of a single tight or diffuse function per type of shell.
# # The tight function is generated from two largest main orbital exponents of the same shell.
# # The diffuse function is generated from two smallest main orbital exponents of the same shell.
# ###############################################################################################
def process_tight_diffuse_func(orbital_exponents_np, geom_mean, orbital_type):
    if orbital_exponents_np.size not in {0, 1}:
        tight_exponent = ((orbital_exponents_np[1:] ** 2) / orbital_exponents_np[:-1])
        exponents_summary_np = concatenate_and_short(geom_mean.flat, tight_exponent[-1,].flat, kind='mergesort')
        diffuse_exponent = ((orbital_exponents_np[:-1] ** 2) / orbital_exponents_np[1:])
        exponents_summary_np = concatenate_and_short(exponents_summary_np.flat, diffuse_exponent[0].flat, kind='mergesort')
        print(f'\n {orbital_type} orbital exponents, tight and a diffuse exponent:\n', exponents_summary_np)
    else:
        exponents_summary_np = geom_mean.flatten()
    return exponents_summary_np

S_exponents_summary_np = process_tight_diffuse_func(S_orbital_exponents_np, geom_mean_S, 'S')

if P_orbital_exponents_np.size not in {0,1}:
    if element_ID not in {'H','He','Li','Be','Na','Mg','K','Ca','Rb','Sr'}:
        tight_exponent_P_np = ((P_orbital_exponents_np[1:] ** 2) / P_orbital_exponents_np[:-1])
        tight_exponent_P_np_2nd = tight_exponent_P_np*4.0
        tight_exponent_P_np_3rd = tight_exponent_P_np*16.0
        tight_exponent_P_np_with_2nd = concatenate_and_short(tight_exponent_P_np_2nd[-1].flat,tight_exponent_P_np[-1,].flat, kind='mergesort')
        tight_exponent_P_np_with_2nd_and_3rd = concatenate_and_short(tight_exponent_P_np_with_2nd.flat,tight_exponent_P_np_3rd[-1,].flat,kind='mergesort')
        P_exponents_summary_np = concatenate_and_short(geom_mean_P.flat,tight_exponent_P_np_with_2nd_and_3rd.flat, kind='mergesort')
        diffuse_exponent_P_np = ((P_orbital_exponents_np[:-1] ** 2) / P_orbital_exponents_np[1:])
        P_exponents_summary_np = concatenate_and_short(P_exponents_summary_np.flat,diffuse_exponent_P_np[0].flat, kind='mergesort')
        print('\n P orbital exponents, tight (+2p) and a diffuse exponent:\n', P_exponents_summary_np)
    else:
        tight_exponent_P_np = ((P_orbital_exponents_np[1:] ** 2) / P_orbital_exponents_np[:-1])
        P_exponents_summary_np = concatenate_and_short(geom_mean_P.flat,tight_exponent_P_np[-1,].flat, kind='mergesort')
        diffuse_exponent_P_np = ((P_orbital_exponents_np[:-1] ** 2) / P_orbital_exponents_np[1:])
        P_exponents_summary_np = concatenate_and_short(P_exponents_summary_np.flat,diffuse_exponent_P_np[0].flat, kind='mergesort')
        print('\n P orbital exponents, tight and a diffuse exponent:\n', P_exponents_summary_np)
else:
    P_exponents_summary_np = geom_mean_P.flatten()

D_exponents_summary_np = process_tight_diffuse_func(D_orbital_exponents_np, geom_mean_D, 'D')
F_exponents_summary_np = process_tight_diffuse_func(F_orbital_exponents_np, geom_mean_F, 'F')
G_exponents_summary_np = process_tight_diffuse_func(G_orbital_exponents_np, geom_mean_G, 'G')
H_exponents_summary_np = process_tight_diffuse_func(H_orbital_exponents_np, geom_mean_H, 'H')
I_exponents_summary_np = process_tight_diffuse_func(I_orbital_exponents_np, geom_mean_I, 'I')
 #############################################################################
# # Addition of a single layer of L-shell exponents by taking the geometric
# mean of the L-1 shell generated exponents. For example, consider the atom
# hydrogen, the AutoCABS with exponents generated by the geometric mean
# plus the tight and diffuse exponent is [4s,3p]. Once an additional layer of
# exponents is added, the AutoCABS is extended to [4s,3p,2d].
# #############################################################################

while(True):
    if D_orbital_exponents_np.size == 0 and P_orbital_exponents_np.size not in {0,1}:
        D_exponents_summary_np = (P_exponents_summary_np[:-1] * P_exponents_summary_np[1:]) ** (1 / 2)
        D_exponents_summary_np = unique_rows(np.vstack(D_exponents_summary_np)).flatten()
        print('Generated D orbital exponents (extended version)\n', D_exponents_summary_np)
        if F_exponents_summary_np.size == 0 and D_exponents_summary_np.size not in {0,1}:
            F_exponents_summary_np = (D_exponents_summary_np[:-1] * D_exponents_summary_np[1:]) ** (1 / 2)
            F_exponents_summary_np = unique_rows(np.vstack(F_exponents_summary_np)).flatten()
            print('Generated F orbital exponents (extended version)\n', F_exponents_summary_np)
            break
        else:
            break
    elif F_exponents_summary_np.size == 0 and D_exponents_summary_np.size not in {0,1}:
        F_exponents_summary_np = (D_exponents_summary_np[:-1] * D_exponents_summary_np[1:]) ** (1 / 2)
        F_exponents_summary_np = unique_rows(np.vstack(F_exponents_summary_np)).flatten()
        print('Generated F orbital exponents (extended version)\n', F_exponents_summary_np)
        if G_exponents_summary_np.size == 0 and F_exponents_summary_np.size not in {0,1}:
            G_exponents_summary_np = (F_exponents_summary_np[:-1] * F_exponents_summary_np[1:]) ** (1 / 2)
            G_exponents_summary_np = unique_rows(np.vstack(G_exponents_summary_np)).flatten()
            print('Generated G orbital exponents (extended version)\n', G_exponents_summary_np)
            break
        else:
            break
    elif G_exponents_summary_np.size == 0 and F_exponents_summary_np.size not in {0,1}:
        G_exponents_summary_np = (F_exponents_summary_np[:-1] * F_exponents_summary_np[1:]) ** (1 / 2)
        G_exponents_summary_np = unique_rows(np.vstack(G_exponents_summary_np)).flatten()
        print('Generated G orbital exponents (extended version)\n', G_exponents_summary_np)
        if H_orbital_exponents_np.size == 0 and G_exponents_summary_np.size not in {0, 1}:
            H_exponents_summary_np = (G_exponents_summary_np[:-1] * G_exponents_summary_np[1:]) ** (1 / 2)
            H_exponents_summary_np = unique_rows(np.vstack(H_exponents_summary_np)).flatten()
            print('Generated H orbital exponents (extended version)\n', H_exponents_summary_np)
            break
        else:
            break
    elif H_orbital_exponents_np.size == 0 and G_exponents_summary_np.size not in {0, 1}:
        H_exponents_summary_np = (G_exponents_summary_np[:-1] * G_exponents_summary_np[1:]) ** (1 / 2)
        H_exponents_summary_np = unique_rows(np.vstack(H_exponents_summary_np)).flatten()
        print('Generated H orbital exponents (extended version)\n', H_exponents_summary_np)
        break
    elif I_orbital_exponents_np.size == 0 and H_exponents_summary_np.size not in {0, 1}:
        I_exponents_summary_np = (H_exponents_summary_np[:-1] * H_exponents_summary_np[1:]) ** (1 / 2)
        I_exponents_summary_np = unique_rows(np.vstack(I_exponents_summary_np)).flatten()
        print('Generated I orbital exponents (extended version)\n', I_exponents_summary_np)
        break
    else:
        break


# Export the AutoCABS to an ORCA-compatible file
write_exponents_summary_to_gbs_file(S_exponents_summary_np,'S_exponents_summary_np','S')
write_exponents_summary_to_gbs_file(P_exponents_summary_np,'P_exponents_summary_np','P')
write_exponents_summary_to_gbs_file(D_exponents_summary_np,'D_exponents_summary_np','D')
write_exponents_summary_to_gbs_file(F_exponents_summary_np,'F_exponents_summary_np','F')
write_exponents_summary_to_gbs_file(G_exponents_summary_np,'G_exponents_summary_np','G')
write_exponents_summary_to_gbs_file(H_exponents_summary_np,'H_exponents_summary_np','H')
write_exponents_summary_to_gbs_file(I_exponents_summary_np,'I_exponents_summary_np','I')   

generated_exponents_summary_gbs = ["S_exponents_summary_np_gbs","P_exponents_summary_np_gbs","D_exponents_summary_np_gbs",
                               "F_exponents_summary_np_gbs","G_exponents_summary_np_gbs","H_exponents_summary_np_gbs","I_exponents_summary_np_gbs"]

with open("%s-AutoCABS.gbs" % element_ID,"w") as outfile:
    outfile.write('%s 0\n' % element_ID)
    for filename in generated_exponents_summary_gbs:
        if os.path.exists(filename) == True:
            with open(filename) as infile:
                contents = infile.read()
                outfile.write(contents)
        if os.path.exists(f'{filename}') == True:
            os.remove(filename)
    outfile.write('****\n')

write_exponents_summary_to_orca_file(S_exponents_summary_np,'S_exponents_summary_np','S')
write_exponents_summary_to_orca_file(P_exponents_summary_np,'P_exponents_summary_np','P')
write_exponents_summary_to_orca_file(D_exponents_summary_np,'D_exponents_summary_np','D')
write_exponents_summary_to_orca_file(F_exponents_summary_np,'F_exponents_summary_np','F')
write_exponents_summary_to_orca_file(G_exponents_summary_np,'G_exponents_summary_np','G')
write_exponents_summary_to_orca_file(H_exponents_summary_np,'H_exponents_summary_np','H')
write_exponents_summary_to_orca_file(I_exponents_summary_np,'I_exponents_summary_np','I') 
generated_exponents_summary_orca = ["S_exponents_summary_np_gbs","P_exponents_summary_np_gbs","D_exponents_summary_np_gbs",
                               "F_exponents_summary_np_gbs","G_exponents_summary_np_gbs","H_exponents_summary_np_gbs","I_exponents_summary_np_gbs"]

# Next, save the generated AutoCABS in an ORCA format
with open("%s-AutoCABS.orca" % element_ID,"w") as outfile:
    outfile.write('$DATA\n\n')
    outfile.write('%s\n' % element_ID)
    for filename in generated_exponents_summary_orca:
        if os.path.exists(filename) == True:
            with open(filename) as infile:
                contents = infile.read()
                outfile.write(contents)
        if os.path.exists(f'{filename}') == True:
            os.remove(filename)
    outfile.write('$END')

# Save the generated AutoCABS in a Turbomole format
write_exponents_summary_to_tm_file(S_exponents_summary_np,'S_exponents_summary_np','S')
write_exponents_summary_to_tm_file(P_exponents_summary_np,'P_exponents_summary_np','P')
write_exponents_summary_to_tm_file(D_exponents_summary_np,'D_exponents_summary_np','D')
write_exponents_summary_to_tm_file(F_exponents_summary_np,'F_exponents_summary_np','F')
write_exponents_summary_to_tm_file(G_exponents_summary_np,'G_exponents_summary_np','G')
write_exponents_summary_to_tm_file(H_exponents_summary_np,'H_exponents_summary_np','H')
write_exponents_summary_to_tm_file(I_exponents_summary_np,'I_exponents_summary_np','I')
generated_exponents_summary_tm = ["S_exponents_summary_np_tm","P_exponents_summary_np_tm","D_exponents_summary_np_tm",
                               "F_exponents_summary_np_tm","G_exponents_summary_np_tm","H_exponents_summary_np_tm","I_exponents_summary_np_tm"]

with open("%s-AutoCABS.0.tm" % element_ID,"w") as outfile:
    outfile.write('$basis\n*\n')
    outfile.write('%s\n#\n*\n' % element_ID)
    for filename in generated_exponents_summary_tm:
        if os.path.exists(filename) == True:
            with open(filename) as infile:
                contents = infile.read()
                outfile.write(contents)
        if os.path.exists(f'{filename}') == True:
            os.remove(filename)
    outfile.write('*\n$end\n')

# Save the generated AutoCABS in a Psi4 format
with open("%s-AutoCABS.gbs" % element_ID) as file1:
    with open("%s-AutoCABS.psi4" % element_ID, "w") as file2:
        file2.write('****\n')
        for line in file1:
            file2.write(line)

# Save the generated AutoCABS in a MOLPRO format
bse.convert_formatted_basis_file("%s-AutoCABS.gbs" % element_ID, "%s-AutoCABS.mpro" % element_ID)

discardthis = "%s-AutoCABS.gbs" % element_ID
if os.path.exists(f'{discardthis}') == True:
    os.remove(discardthis)
discardthisone = '%s.gbs-formatted_basis' % element_ID
if os.path.exists(f'{discardthisone}') == True:
    os.remove(discardthisone)
discardthisone2 = '%s.orca-formatted_basis' % element_ID
if os.path.exists(f'{discardthisone2}') == True:
    os.remove(discardthisone2)
discardthisone3 = '%s.gbs' % element_ID
if os.path.exists(f'{discardthisone3}') == True:
    os.remove(discardthisone3)
print('\nThe AutoCABS was saved in ORCA, Turbomole, Psi4, and MOLPRO compatible formats.')
print('Normal Termination\n')
