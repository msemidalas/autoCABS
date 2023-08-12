
This Python implementation provides an automated approach to generate auxiliary-complementary basis sets (autoCABS) for explicitly correlated F12 calculations.

For a double-zeta orbital basis (cardinal number is D), the script creates an autoCABS by generating geometrically averaged exponents, diffuse and tight functions plus two tight p exponents for the p-block elements, and two additional layers of exponents. For larger orbital basis sets, n = T, Q, 5, the program incorporates the geometrically averaged exponents, the diffuse functions, tight functions without additional p exponents, and only a single layer of exponents.

To use the program, follow these steps:

- Prepare a basis set file in ORCA or MOLPRO format for each element, naming them as 'Element.orca' or 'Element.mpro', respectively (e.g., 'c.orca' or 'c.mpro' for carbon).

- Place three files in a single folder: the autoCABS generation script, 'my_functions.py' with the function definitions, and the orbital basis set file (in either ORCA or MOLPRO format).

- To generate an autoCABS for n = D, execute the following command:

```python3 AutoCABS-generation-geom_mean_tight_extended_2p_and_diffuse_plus_two_layers.py```

To generate an autoCABS from a larger orbital basis set with a cardinal number n larger than D, use

```python3 AutoCABS-generation-geom_mean_tight_and_diffuse_plus_single_layer.py```

The autoCABS is exported in four different formats, MOLPRO, ORCA, TURBOMOLE, and Psi4.

To generate autoCABS for multiple elements, repeat the above steps for each element. A simple bash script is provided to assist in combining the generated autoCABS files into a single file.

Feel free to contribute, provide feedback, or report issues.

Emmanouil Semidalas
Email: msemidalas@yahoo.com

**References** 

Semidalas E, Martin JML. Automatic generation of complementary auxiliary basis sets for explicitly correlated methods. J Comput Chem. 2022;43(25):1690-1700. doi:10.1002/jcc.26970
