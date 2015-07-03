Black phosphorus Top-Of-Barrier code with EK from vasp HSE calculation using Wannier interpolation

change pORn='n' or pORn='p' in input_file.m to choose n-type or p-type device (using EC or EV)
change n_Angl number for different transport angle (current setting for x (0) and y (90) transport directions only)

Calculation in Vg done from 0 to 0.8V and the results are interpolated in plot_IV2.m to align all I-OFF to 0.1uA/um

Figures for paper plotted using plot_IV2.m