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

# Reads and edits a MOLPRO basis (e.g.,c.mpro) and converts it to an ORCA-compatible one for further operations.
for f in glob.glob("-AutoCABS*":
    os.remove(f)

molpro_filenames = glob.glob('**/*.mpro',recursive=True)
filename_molpro_original = len(molpro_filenames)
if filename_molpro_original is not 0:
    print('Basis set file in a MOLPRO format (.mpro) was found.')
element_ID = (molpro_filenames[0]).split(".")[0].title()
print('Chosen element is: %s' % element_ID)
formatted_basis_mpro = BasisSet.from_file('%s' % molpro_filenames[0], fmt='molpro', name='user-def-basis')
formatted_basis_gbs = formatted_basis_mpro.to_gaussian()
with open('%s.gbs' % element_ID,"w") as outfile:
    outfile.write(formatted_basis_gbs)
read_and_format_gbs_basis_set('%s.gbs' % element_ID)
formatted_basis = '%s.gbs-formatted_basis' % element_ID

#orca_filenames = glob.glob('**/*.orca',recursive=True)
#filename_original = len(orca_filenames)
#if filename_original == 0:
#    print('Basis set file was not found.')
#    sys.exit()
#print((orca_filenames[0]).split(".")[0].title())
#element_ID = (orca_filenames[0]).split(".")[0].title()
#print('Chosen element is:%s' % element_ID)
#formatted_basis = '%s-formatted_basis' % orca_filenames[0]
#read_and_format_orca_basis_set(orca_filenames[0])

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
df['Contraction'] = np.where(df['Contraction_Coeff'] == 1,1,0)
pd.set_option('display.max_columns',None)

# Group exponents in dataframes by shells.
S_orbital_exponents = df[df['S_Type_Exp']==1]
P_orbital_exponents = df[df['P_Type_Exp']==1]
D_orbital_exponents = df[df['D_Type_Exp']==1]
F_orbital_exponents = df[df['F_Type_Exp']==1]
G_orbital_exponents = df[df['G_Type_Exp']==1]
H_orbital_exponents = df[df['H_Type_Exp']==1]

# Group exponents based on whether their contraction coeffs. are equal to one or not.
# Contraction coeff. equals to 1
S_orbital_exponents_with_contraction_coeff_1 = S_orbital_exponents[S_orbital_exponents['Contraction']==1]
S_orbital_exponents_with_contraction_coeff_1_np = np.array(S_orbital_exponents_with_contraction_coeff_1['Orbital_Exponent'],dtype="float64")
if S_orbital_exponents_with_contraction_coeff_1_np.size != 0:
    S_orbital_exponents_with_contraction_coeff_1_np = np.vstack(S_orbital_exponents_with_contraction_coeff_1_np)
    S_orbital_exponents_with_contraction_coeff_1_np = unique_rows(S_orbital_exponents_with_contraction_coeff_1_np)
    S_orbital_exponents_with_contraction_coeff_1_np = np.sort(S_orbital_exponents_with_contraction_coeff_1_np, axis=-1)
#Contraction coeff. does not equal to 1
S_orbital_exponents_with_contraction_coeff_not_1 = S_orbital_exponents[S_orbital_exponents['Contraction']!=1]
S_orbital_exponents_with_contraction_coeff_not_1_np = np.array(S_orbital_exponents_with_contraction_coeff_not_1['Orbital_Exponent'],dtype="float64")
if S_orbital_exponents_with_contraction_coeff_not_1_np.size != 0:
    S_orbital_exponents_with_contraction_coeff_not_1_np = np.vstack(S_orbital_exponents_with_contraction_coeff_not_1_np)
    S_orbital_exponents_with_contraction_coeff_not_1_np = unique_rows(S_orbital_exponents_with_contraction_coeff_not_1_np)
    S_orbital_exponents_with_contraction_coeff_not_1_np = np.sort(S_orbital_exponents_with_contraction_coeff_not_1_np, axis=-1)
print('\n S main orbital exponents with contraction coefficient equal to 1:\n',S_orbital_exponents_with_contraction_coeff_1_np)
print('\n S main orbital exponents with contraction coefficient not equal to 1:\n',S_orbital_exponents_with_contraction_coeff_not_1_np)

#P shells
#Contraction coeff. equals to 1
P_orbital_exponents_with_contraction_coeff_1 = P_orbital_exponents[P_orbital_exponents['Contraction']==1]
P_orbital_exponents_with_contraction_coeff_1_np = np.array(P_orbital_exponents_with_contraction_coeff_1['Orbital_Exponent'],dtype="float64")
if P_orbital_exponents_with_contraction_coeff_1_np.size != 0:
    P_orbital_exponents_with_contraction_coeff_1_np = np.vstack(P_orbital_exponents_with_contraction_coeff_1_np)
    P_orbital_exponents_with_contraction_coeff_1_np = unique_rows(P_orbital_exponents_with_contraction_coeff_1_np)
    P_orbital_exponents_with_contraction_coeff_1_np = np.sort(P_orbital_exponents_with_contraction_coeff_1_np, axis=-1)
#Contraction coeff. does not equal to 1
P_orbital_exponents_with_contraction_coeff_not_1 = P_orbital_exponents[P_orbital_exponents['Contraction']!=1]
P_orbital_exponents_with_contraction_coeff_not_1_np = np.array(P_orbital_exponents_with_contraction_coeff_not_1['Orbital_Exponent'],dtype="float64")
if P_orbital_exponents_with_contraction_coeff_not_1_np.size != 0:
    P_orbital_exponents_with_contraction_coeff_not_1_np = np.vstack(P_orbital_exponents_with_contraction_coeff_not_1_np)
    P_orbital_exponents_with_contraction_coeff_not_1_np = unique_rows(P_orbital_exponents_with_contraction_coeff_not_1_np)
    P_orbital_exponents_with_contraction_coeff_not_1_np = np.sort(P_orbital_exponents_with_contraction_coeff_not_1_np, axis=-1)
print('\n P main orbital exponents with contraction coefficient equal to 1:\n',P_orbital_exponents_with_contraction_coeff_1_np)
print('\n P main orbital exponents with contraction coefficient not equal to 1:\n',P_orbital_exponents_with_contraction_coeff_not_1_np)

# D shells
#Contraction coeff. equals to 1
D_orbital_exponents_with_contraction_coeff_1 = D_orbital_exponents[D_orbital_exponents['Contraction']==1]
D_orbital_exponents_with_contraction_coeff_1_np = np.array(D_orbital_exponents_with_contraction_coeff_1['Orbital_Exponent'],dtype="float64")
if D_orbital_exponents_with_contraction_coeff_1_np.size != 0:
    D_orbital_exponents_with_contraction_coeff_1_np = np.vstack(D_orbital_exponents_with_contraction_coeff_1_np)
    D_orbital_exponents_with_contraction_coeff_1_np = unique_rows(D_orbital_exponents_with_contraction_coeff_1_np)
    D_orbital_exponents_with_contraction_coeff_1_np = np.sort(D_orbital_exponents_with_contraction_coeff_1_np, axis=-1)
#Contraction coeff. does not equal to 1
D_orbital_exponents_with_contraction_coeff_not_1 = D_orbital_exponents[D_orbital_exponents['Contraction']!=1]
D_orbital_exponents_with_contraction_coeff_not_1_np = np.array(D_orbital_exponents_with_contraction_coeff_not_1['Orbital_Exponent'],dtype="float64")
if D_orbital_exponents_with_contraction_coeff_not_1_np.size != 0:
    D_orbital_exponents_with_contraction_coeff_not_1_np = np.vstack(D_orbital_exponents_with_contraction_coeff_not_1_np)
    D_orbital_exponents_with_contraction_coeff_not_1_np = unique_rows(D_orbital_exponents_with_contraction_coeff_not_1_np)
    D_orbital_exponents_with_contraction_coeff_not_1_np = np.sort(D_orbital_exponents_with_contraction_coeff_not_1_np, axis=-1)
print('\n D main orbital exponents with contraction coefficient equal to 1:\n',D_orbital_exponents_with_contraction_coeff_1_np)
print('\n D main orbital exponents with contraction coefficient not equal to 1:\n',D_orbital_exponents_with_contraction_coeff_not_1_np)

# F shells
# Contraction coeff. equals to 1
F_orbital_exponents_with_contraction_coeff_1 = F_orbital_exponents[F_orbital_exponents['Contraction']==1]
F_orbital_exponents_with_contraction_coeff_1_np = np.array(F_orbital_exponents_with_contraction_coeff_1['Orbital_Exponent'],dtype="float64")
if F_orbital_exponents_with_contraction_coeff_1_np.size != 0:
    F_orbital_exponents_with_contraction_coeff_1_np = np.vstack(F_orbital_exponents_with_contraction_coeff_1_np)
    F_orbital_exponents_with_contraction_coeff_1_np = unique_rows(F_orbital_exponents_with_contraction_coeff_1_np)
    F_orbital_exponents_with_contraction_coeff_1_np = np.sort(F_orbital_exponents_with_contraction_coeff_1_np, axis=-1)
# Contraction coeff. does not equal to 1
F_orbital_exponents_with_contraction_coeff_not_1 = F_orbital_exponents[F_orbital_exponents['Contraction']!=1]
F_orbital_exponents_with_contraction_coeff_not_1_np = np.array(F_orbital_exponents_with_contraction_coeff_not_1['Orbital_Exponent'],dtype="float64")
if F_orbital_exponents_with_contraction_coeff_not_1_np.size != 0:
    F_orbital_exponents_with_contraction_coeff_not_1_np = np.vstack(F_orbital_exponents_with_contraction_coeff_not_1_np)
    F_orbital_exponents_with_contraction_coeff_not_1_np = unique_rows(F_orbital_exponents_with_contraction_coeff_not_1_np)
    F_orbital_exponents_with_contraction_coeff_not_1_np = np.sort(F_orbital_exponents_with_contraction_coeff_not_1_np, axis=-1)
print('\n F main orbital exponents with contraction coefficient equal to 1:\n',F_orbital_exponents_with_contraction_coeff_1_np)
print('\n F main orbital exponents with contraction coefficient not equal to 1:\n',F_orbital_exponents_with_contraction_coeff_not_1_np)

# G shells
G_orbital_exponents_with_contraction_coeff_1 = G_orbital_exponents[G_orbital_exponents['Contraction']==1]
G_orbital_exponents_with_contraction_coeff_1_np = np.array(G_orbital_exponents_with_contraction_coeff_1['Orbital_Exponent'],dtype="float64")
if G_orbital_exponents_with_contraction_coeff_1_np.size != 0:
    G_orbital_exponents_with_contraction_coeff_1_np = np.vstack(G_orbital_exponents_with_contraction_coeff_1_np)
    G_orbital_exponents_with_contraction_coeff_1_np = unique_rows(G_orbital_exponents_with_contraction_coeff_1_np)
    G_orbital_exponents_with_contraction_coeff_1_np = np.sort(G_orbital_exponents_with_contraction_coeff_1_np, axis=-1)
#Contraction coeff. does not equal to 1
G_orbital_exponents_with_contraction_coeff_not_1 = G_orbital_exponents[G_orbital_exponents['Contraction']!=1]
G_orbital_exponents_with_contraction_coeff_not_1_np = np.array(G_orbital_exponents_with_contraction_coeff_not_1['Orbital_Exponent'],dtype="float64")
if G_orbital_exponents_with_contraction_coeff_not_1_np.size != 0:
    G_orbital_exponents_with_contraction_coeff_not_1_np = np.vstack(G_orbital_exponents_with_contraction_coeff_not_1_np)
    G_orbital_exponents_with_contraction_coeff_not_1_np = unique_rows(G_orbital_exponents_with_contraction_coeff_not_1_np)
    G_orbital_exponents_with_contraction_coeff_not_1_np = np.sort(G_orbital_exponents_with_contraction_coeff_not_1_np, axis=-1)
print('\n G main orbital exponents with contraction coefficient equal to 1:\n',G_orbital_exponents_with_contraction_coeff_1_np)
print('\n G main orbital exponents with contraction coefficient not equal to 1:\n',G_orbital_exponents_with_contraction_coeff_not_1_np)

#H shells
H_orbital_exponents_with_contraction_coeff_1 = H_orbital_exponents[H_orbital_exponents['Contraction']==1]
H_orbital_exponents_with_contraction_coeff_1_np = np.array(H_orbital_exponents_with_contraction_coeff_1['Orbital_Exponent'],dtype="float64")
if H_orbital_exponents_with_contraction_coeff_1_np.size != 0:
    H_orbital_exponents_with_contraction_coeff_1_np = np.vstack(H_orbital_exponents_with_contraction_coeff_1_np)
    H_orbital_exponents_with_contraction_coeff_1_np = unique_rows(H_orbital_exponents_with_contraction_coeff_1_np)
    H_orbital_exponents_with_contraction_coeff_1_np = np.sort(H_orbital_exponents_with_contraction_coeff_1_np, axis=-1)
#Contraction coeff. does not equal to 1
H_orbital_exponents_with_contraction_coeff_not_1 = H_orbital_exponents[H_orbital_exponents['Contraction']!=1]
H_orbital_exponents_with_contraction_coeff_not_1_np = np.array(H_orbital_exponents_with_contraction_coeff_not_1['Orbital_Exponent'],dtype="float64")
if H_orbital_exponents_with_contraction_coeff_not_1_np.size != 0:
    H_orbital_exponents_with_contraction_coeff_not_1_np = np.vstack(H_orbital_exponents_with_contraction_coeff_not_1_np)
    H_orbital_exponents_with_contraction_coeff_not_1_np = unique_rows(H_orbital_exponents_with_contraction_coeff_not_1_np)
    H_orbital_exponents_with_contraction_coeff_not_1_np = np.sort(H_orbital_exponents_with_contraction_coeff_not_1_np, axis=-1)
print('\n H main orbital exponents with contraction coefficient equal to 1:\n',H_orbital_exponents_with_contraction_coeff_1_np)
print('\n H main orbital exponents with contraction coefficient not equal to 1:\n',H_orbital_exponents_with_contraction_coeff_not_1_np)

#############################################################################################################################
##### Grab exponents with a contraction coeff. equal to 1 plus an outermost exponent with contraction coeff not equal to 1 ##
#############################################################################################################################
# S exponents
# Find outermost exponent with contraction coefficient not equal to 1
S_orbital_exponents_outermost_decontracted = np.setdiff1d(S_orbital_exponents_with_contraction_coeff_not_1_np.flatten(),
                                                          S_orbital_exponents_with_contraction_coeff_1_np.flatten())
if S_orbital_exponents_outermost_decontracted.size != 0:
    S_orbital_exponents_np = concatenate_and_short(S_orbital_exponents_with_contraction_coeff_1_np.flat,
                                                   S_orbital_exponents_outermost_decontracted[0,].flat,kind='mergesort')
    print('\n S orbital exponents and a single outermost decontracted exponent\n',S_orbital_exponents_np)
else:
    S_orbital_exponents_np = S_orbital_exponents_with_contraction_coeff_1_np
    print('\n S orbital exponents are all with contraction coeff. equal to 1\n',S_orbital_exponents_np)

# P exponents
# Find outermost exponent with contraction coefficient not equal to 1
P_orbital_exponents_outermost_decontracted = np.setdiff1d(P_orbital_exponents_with_contraction_coeff_not_1_np.flatten(),
                                                          P_orbital_exponents_with_contraction_coeff_1_np.flatten())
if P_orbital_exponents_outermost_decontracted.size != 0:
    P_orbital_exponents_np = concatenate_and_short(P_orbital_exponents_with_contraction_coeff_1_np.flat,
                                                   P_orbital_exponents_outermost_decontracted[0,].flat,kind='mergesort')
    print('\n P orbital exponents and a single outermost decontracted exponent\n',P_orbital_exponents_np)
else:
    P_orbital_exponents_np = P_orbital_exponents_with_contraction_coeff_1_np
    print('\n P orbital exponents are all with contraction coeff. equal to 1\n',P_orbital_exponents_np)

# D exponents
# Find outermost exponent with contraction coefficient not equal to 1
D_orbital_exponents_outermost_decontracted = np.setdiff1d(D_orbital_exponents_with_contraction_coeff_not_1_np.flatten(),
                                                          D_orbital_exponents_with_contraction_coeff_1_np.flatten())
if D_orbital_exponents_outermost_decontracted.size != 0:
    D_orbital_exponents_np = concatenate_and_short(D_orbital_exponents_with_contraction_coeff_1_np.flat,
                                                   D_orbital_exponents_outermost_decontracted[0,].flat,kind='mergesort')
    print('\n D orbital exponents and a single outermost decontracted exponent\n',D_orbital_exponents_np)
else:
    D_orbital_exponents_np = D_orbital_exponents_with_contraction_coeff_1_np
    print('\n D orbital exponents are all with contraction coeff. equal to 1\n',D_orbital_exponents_np)

# F exponents
# Find outermost exponent with contraction coefficient not equal to 1
F_orbital_exponents_outermost_decontracted = np.setdiff1d(F_orbital_exponents_with_contraction_coeff_not_1_np.flatten(),
                                                          F_orbital_exponents_with_contraction_coeff_1_np.flatten())
if F_orbital_exponents_outermost_decontracted.size != 0:
    F_orbital_exponents_np = concatenate_and_short(F_orbital_exponents_with_contraction_coeff_1_np.flat,
                                                   F_orbital_exponents_outermost_decontracted[0,].flat,kind='mergesort')
    print('\n F orbital exponents and a single outermost decontracted exponent\n',F_orbital_exponents_np)
else:
    F_orbital_exponents_np = F_orbital_exponents_with_contraction_coeff_1_np
    print('\n F orbital exponents are all with contraction coeff. equal to 1\n',F_orbital_exponents_np)

# G exponents
# Find outermost exponent with contraction coefficient not equal to 1
G_orbital_exponents_outermost_decontracted = np.setdiff1d(G_orbital_exponents_with_contraction_coeff_not_1_np.flatten(),
                                                          G_orbital_exponents_with_contraction_coeff_1_np.flatten())
if G_orbital_exponents_outermost_decontracted.size != 0:
    G_orbital_exponents_np = concatenate_and_short(G_orbital_exponents_with_contraction_coeff_1_np.flat,
                                                   G_orbital_exponents_outermost_decontracted[0,].flat,kind='mergesort')
    print('\n G orbital exponents and a single outermost decontracted exponent\n',G_orbital_exponents_np)
else:
    G_orbital_exponents_np = G_orbital_exponents_with_contraction_coeff_1_np
    print('\n G orbital exponents are all with contraction coeff. equal to 1\n',G_orbital_exponents_np)

# H exponents
# Find outermost exponent with contraction coefficient not equal to 1
H_orbital_exponents_outermost_decontracted = np.setdiff1d(H_orbital_exponents_with_contraction_coeff_not_1_np.flatten(),
                                                          H_orbital_exponents_with_contraction_coeff_1_np.flatten())
if H_orbital_exponents_outermost_decontracted.size != 0:
    H_orbital_exponents_np = concatenate_and_short(H_orbital_exponents_with_contraction_coeff_1_np.flat,
                                                   H_orbital_exponents_outermost_decontracted[0,].flat,kind='mergesort')
    print('\n H orbital exponents and a single outermost decontracted exponent\n',H_orbital_exponents_np)
else:
    H_orbital_exponents_np = H_orbital_exponents_with_contraction_coeff_1_np
    print('\n H orbital exponents are all with contraction coeff. equal to 1\n',H_orbital_exponents_np)

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

# Summary with generated exponents
print('\nSummary of generated AutoCABS exponents')
print('S\n', np.sort(geom_mean_S)[::-1])
print('P\n', np.sort(geom_mean_P)[::-1])
print('D\n', np.sort(geom_mean_D)[::-1])
print('F\n', np.sort(geom_mean_F)[::-1])
print('G\n', np.sort(geom_mean_G)[::-1])
print('H\n', np.sort(geom_mean_H)[::-1])

# ##############################################################################################
# # Addition of a single tight or diffuse function per type of shell.
# # The tight function is generated from two largest main orbital exponents of the same shell.
# # The diffuse function is generated from two smallest main orbital exponents of the same shell.
# ###############################################################################################
if S_orbital_exponents_np.size not in {0,1}:
    tight_exponent_S_np = ((S_orbital_exponents_np[1:] ** 2) / S_orbital_exponents_np[:-1])
    S_exponents_summary_np = concatenate_and_short(geom_mean_S.flat,tight_exponent_S_np[-1,].flat, kind='mergesort')
    diffuse_exponent_S_np = ((S_orbital_exponents_np[:-1] ** 2) / S_orbital_exponents_np[1:])
    S_exponents_summary_np = concatenate_and_short(S_exponents_summary_np.flat,diffuse_exponent_S_np[0].flat, kind='mergesort')
    print('\n S orbital exponents, tight and a diffuse exponent:\n', S_exponents_summary_np)
else:
    S_exponents_summary_np = geom_mean_S.flatten()

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

if D_orbital_exponents_np.size not in {0,1}:
    tight_exponent_D_np = ((D_orbital_exponents_np[1:] ** 2) / D_orbital_exponents_np[:-1])
    D_exponents_summary_np = concatenate_and_short(geom_mean_D.flat,tight_exponent_D_np[-1,].flat, kind='mergesort')
    diffuse_exponent_D_np = ((D_orbital_exponents_np[:-1] ** 2) / D_orbital_exponents_np[1:])
    D_exponents_summary_np = concatenate_and_short(D_exponents_summary_np.flat,diffuse_exponent_D_np[0].flat, kind='mergesort')
    print('\n D orbital exponents, tight and a diffuse exponent:\n', D_exponents_summary_np)
else:
    D_exponents_summary_np = geom_mean_D.flatten()

if F_orbital_exponents_np.size not in {0,1}:
    tight_exponent_F_np = ((F_orbital_exponents_np[1:] ** 2) / F_orbital_exponents_np[:-1])
    F_exponents_summary_np = concatenate_and_short(geom_mean_F.flat,tight_exponent_F_np[-1,].flat, kind='mergesort')
    diffuse_exponent_F_np = ((F_orbital_exponents_np[:-1] ** 2) / F_orbital_exponents_np[1:])
    F_exponents_summary_np = concatenate_and_short(F_exponents_summary_np.flat,diffuse_exponent_F_np[0].flat, kind='mergesort')
    print('\n F orbital exponents, tight and a diffuse exponent:\n', F_exponents_summary_np)
else:
    F_exponents_summary_np = geom_mean_F.flatten()

if G_orbital_exponents_np.size not in {0,1}:
    tight_exponent_G_np = ((G_orbital_exponents_np[1:] ** 2) / G_orbital_exponents_np[:-1])
    G_exponents_summary_np = concatenate_and_short(geom_mean_G.flat,tight_exponent_G_np[-1,].flat, kind='mergesort')
    diffuse_exponent_G_np = ((G_orbital_exponents_np[:-1] ** 2) / G_orbital_exponents_np[1:])
    G_exponents_summary_np = concatenate_and_short(G_exponents_summary_np.flat,diffuse_exponent_G_np[0].flat, kind='mergesort')
    print('\n G orbital exponents, tight and a diffuse exponent:\n', G_exponents_summary_np)
else:
    G_exponents_summary_np = geom_mean_G.flatten()

if H_orbital_exponents_np.size not in {0,1}:
    tight_exponent_H_np = ((H_orbital_exponents_np[1:] ** 2) / H_orbital_exponents_np[:-1])
    H_exponents_summary_np = concatenate_and_short(geom_mean_H.flat,tight_exponent_H_np[-1,].flat, kind='mergesort')
    diffuse_exponent_H_np = ((H_orbital_exponents_np[:-1] ** 2) / H_orbital_exponents_np[1:])
    H_exponents_summary_np = concatenate_and_short(H_exponents_summary_np.flat,diffuse_exponent_H_np[0].flat, kind='mergesort')
    print('\n H orbital exponents, tight and a diffuse exponent:\n', H_exponents_summary_np)
else:
    H_exponents_summary_np = geom_mean_H.flatten()

# ############################################################################
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
    else:
        break

# Export the AutoCABS to an ORCA-compatible file
write_exponents_summary_to_gbs_file(S_exponents_summary_np,'S_exponents_summary_np','S')
write_exponents_summary_to_gbs_file(P_exponents_summary_np,'P_exponents_summary_np','P')
write_exponents_summary_to_gbs_file(D_exponents_summary_np,'D_exponents_summary_np','D')
write_exponents_summary_to_gbs_file(F_exponents_summary_np,'F_exponents_summary_np','F')
write_exponents_summary_to_gbs_file(G_exponents_summary_np,'G_exponents_summary_np','G')
write_exponents_summary_to_gbs_file(H_exponents_summary_np,'H_exponents_summary_np','H')

generated_exponents_summary_gbs = ["S_exponents_summary_np_gbs","P_exponents_summary_np_gbs","D_exponents_summary_np_gbs",
                               "F_exponents_summary_np_gbs","G_exponents_summary_np_gbs","H_exponents_summary_np_gbs"]

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
generated_exponents_summary_orca = ["S_exponents_summary_np_gbs","P_exponents_summary_np_gbs","D_exponents_summary_np_gbs",
                               "F_exponents_summary_np_gbs","G_exponents_summary_np_gbs","H_exponents_summary_np_gbs"]

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

generated_exponents_summary_tm = ["S_exponents_summary_np_tm","P_exponents_summary_np_tm","D_exponents_summary_np_tm",
                               "F_exponents_summary_np_tm","G_exponents_summary_np_tm","H_exponents_summary_np_tm"]

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

os.remove("%s-AutoCABS.gbs" % element_ID)
os.remove('%s.gbs-formatted_basis' % element_ID)
os.remove('%s.gbs' % element_ID)
print('\nThe AutoCABS was saved in ORCA, Turbomole, Psi4, and MOLPRO compatible formats.')
print('Normal Termination\n')
