%num_bands=2;        % Actually valley degeneracy
pORn='p';           % options: 'p','n'
tox = 3;          % Unit: nm
T = 300;
epsox = 25;

%Ef =1.18; % Position of Fermi level sets the off current (or Vt)
alphag = 1;%0.88;
alphad = 0;%0.035;

%Voltage (gate or drain) steps
Vd_i=0.6;
Vd_f=0.6;
n_Vd=1;

Vg_i=0;
Vg_f=0.8;
n_Vg=41;

% Source to drain transport direction in "degree"
Angl_i=0; 
Angl_f=90; 
n_Angl=2;
