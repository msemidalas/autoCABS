
This Python implementation generates an auxiliary-complementary basis set (autoCABS) and is available at Github (see https://github.com/msemidalas/autoCABS). 

The geometrically averaged exponents, the diffuse and tight functions plus two tight p exponents for the p-block elements, and two additional layers of exponents are added to the autoCABS when the cardinal number of the orbital basis set is D. For larger orbital basis sets, n = T, Q, 5, the program similarly adds the geometrically averaged exponents, the diffuse functions, tight functions without additional p exponents, and one single layer of exponents to the generated autoCABS.

For input instructions, you use a basis set file in an ORCA or MOLPRO format for one element at a time, and name it as either ‘Element.orca’ or ‘Element.mpro’ where element refers to the atomic symbol, (for the atom of carbon, e.g., c.orca or c.mpro).
Make sure that you place the three files, the autoCABS generation script, the my_functions.py with the function definitions, and the orbital basis set file (either in an Orca or Molpro format) in a single folder.

Then to generate an autoCABS for n = D, simply type

python3 AutoCABS-generation-geom_mean_tight_extended_2p_and_diffuse_plus_two_layers.py

To generate an autoCABS from a larger orbital basis set with a cardinal number n larger than D,

python3 AutoCABS-generation-geom_mean_tight_and_diffuse_plus_single_layer.py

The autoCABS is exported in four different formats, MOLPRO, ORCA, TURBOMOLE, and Psi4.

Following the procedure described above, the autoCABS can be generated for the other elements and then all of them should be combined into a single file, there is a simple bash script that assists in that.

