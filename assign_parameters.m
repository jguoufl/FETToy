
%epsox = 25;
m0=9.1095e-31;
hbar = 1.05458e-34;
q = 1.60218e-19;
eps0 = 8.85418E-12;
kb = 1.38066e-23;

% thermal voltage in eV
kT = kb*T/q;
% in meter
tox=tox*1e-9;
% the variables must be E_all (2D bandstructure), kx_plot (in /m)

switch pORn
    case 'n'
        %load BP_mono_EC.mat
        
        load BP_mono_dense_EC.mat
        kx_plot=kx_plot*(2*pi/4.5694e-10);
        ky_plot=ky_plot*(2*pi/3.3255e-10);
        %}
        %load MoS2_EC.mat
    case 'p'
        %load BP_mono_EV.mat
        
        load BP_mono_dense_EV.mat
        kx_plot=kx_plot*(2*pi/4.5694e-10);
        ky_plot=ky_plot*(2*pi/3.3255e-10);
        %}
        %load MoS2_EV.mat
end

%{
load TwoD_Ek
ky_plot=kx_plot;
%}

minE=min(min(abs(E_all)));
Ef=minE-0.25;
fprintf('Ef=%.4f\n',Ef);

sizeE_all=size(E_all);
nkx=sizeE_all(1);
dkx=abs(kx_plot(nkx)-kx_plot(nkx-1));
nky=sizeE_all(2);
dky=abs(ky_plot(nky)-ky_plot(nky-1));
n_states=2*(dkx/(2*pi))*(dky/(2*pi));
% IMPORTANT: bands are always assumed spin degenerate, 
%   therefore 2* used here
% nE_all=sizeE_all(3); % number of subbands

num_bands=1;
E=nan*ones(nkx,nkx,num_bands);
ferm=nan*ones(nkx,nkx,num_bands); % Equillibrium fermi distribution
Vel_x=E;
Vel_y=E;
Vel=E;
mu=nan*ones(nkx,nkx,num_bands);

switch pORn
    case 'n'
        E=squeeze(E_all(:,:,1:num_bands));
    case 'p'
        %E=-squeeze(E_all(:,:,(nE_all-num_bands+1:nE_all)));
        E=squeeze(E_all(:,:,1:num_bands));
        % the negative sign ensures that all band look electron like       
end



N0ii=nan*ones(num_bands,1);
tic
for ii=1:num_bands
    [Vel_x(:,:,ii),Vel_y(:,:,ii)]=...
        gradient(q*(1/hbar)*squeeze(E(:,:,ii)),dkx,dky);
    ferm(:,:,ii)=(1./(1+exp((squeeze(E(:,:,ii))-Ef*ones(nkx,nkx))/kT)));
    f_tmp=ferm(:,:,ii);
    A_tmp=isnan(f_tmp);
    A_tmp=find(A_tmp==1);
    f_tmp(A_tmp)=0; 
    N0ii(ii) = sum(sum( (n_states*f_tmp)));
    % note: no spin or valley degeneracy hardwired
end
toc
N0 = sum(N0ii);

% Transport orientation information
%Angl=(pi/180)*(linspace(Angl_i, Angl_f, n_Angl))';
Angl=(linspace(Angl_i, Angl_f, n_Angl))';
%Angl(Angl==0)=1e-3;
%Angl(Angl==90)=90-1e-3;
% Transport orientation w.r.t. [100] direction, unit: radian

% Gate and drain biases
Vd=(linspace(Vd_i, Vd_f, n_Vd))';
if (Vd(1)==0)
    Vd(1)=Vd(1)+1e-3;
end
Vg=(linspace(Vg_i, Vg_f, n_Vg))';

Nll=nan*ones(n_Vd,n_Vg,n_Angl,num_bands);
dN=nan*ones(n_Vd,n_Vg,n_Angl,num_bands);

% the factor of 2 here assumes a double gate device
Cins = 2*epsox*eps0/tox;

% charging energy
C0 = Cins/alphag;
U0 = q/C0;


%Number of electrons in transistor at equilibrium
