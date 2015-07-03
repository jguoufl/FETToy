clear all
close all

style={'b-o', 'b-x', 'r--o', 'r--x', 'b-^', 'r--^'};
%{
files={'BP_x_n.mat';
    'BP_x_p.mat';
    'BP_y_n.mat';
    'BP_y_p.mat'};

figure('units','inches','position',[1 1 4 3]);
axes('units','inches','position',[0.6 0.45 3 2.4]);
for idx=1:4
load(files{idx},'Vd','Vg','I');
vg=Vg;
vd=Vd(end);
curr=I(end,:);

v_off=interp1(log10(curr),vg,-1);
v_on=v_off+vd;
vg_new=linspace(v_off,v_on,21);
curr_new=10.^interp1(vg,log10(curr),vg_new);
vg_new=vg_new-v_off;
semilogy(vg_new,curr_new,style{idx})
hold on;
end
set(gca,'box','off')
xlim([0 vd]);
xlabel('V_G (V)')
ylabel('Current (\muA/\mum)')

axes('units','inches','position',[0.6 0.45 3 2.4]);
for idx=1:4
load(files{idx},'Vd','Vg','I');
vg=Vg;
vd=Vd(end);
curr=I(end,:);

v_off=interp1(log10(curr),vg,-1);
v_on=v_off+vd;
vg_new=linspace(v_off,v_on,21);
curr_new=10.^interp1(vg,log10(curr),vg_new);
vg_new=vg_new-v_off;
plot(vg_new,curr_new,style{idx})
hold on;
end
set(gca,'color','none','yaxislocation','right','xticklabel',[])
xlim([0 vd]);
ylim([0 6000]);
%}

F1=figure('units','inches','position',[1 1 4 3.5]);
axes('units','inches','position',[0.6 2.25 1.2 1.15]);
load data1/BP_n_angle_dense Vg Angl I NQ Uscf_mat Ef E_all
NQ=squeeze(NQ);
Vel_old=I./NQ;
vd=0.6;
vg=Vg;
curr_old=I;
curr_new=zeros(length(Angl),21);
Vel_new_n=zeros(length(Angl),21);
delay_n=cell(length(Angl),1);
Uscf_Vg=zeros(length(Angl),21);
for idx=1:length(Angl)
    curr=curr_old(:,idx);
    v_off=interp1(log10(curr),vg,-1);
    v_on=v_off+vd;
    vg_new=linspace(v_off,v_on,21);
    curr_new(idx,:)=(10.^interp1(vg,log10(curr),vg_new));
    vel=Vel_old(:,idx);
    Vel_new_n(idx,:)=(interp1(vg,vel,vg_new));
    NQ_tmp=NQ(:,idx);
    delay=zeros(11,2);
    for idx_vg=1:11
        dQ=NQ_tmp(idx_vg+[0 30]);
        dI=curr(idx_vg+[0 30]);
        delay(idx_vg,:)=[(dI(2)/dI(1)) (dQ(2)-dQ(1))/(dI(2)-dI(1))];
    end
    delay_n{idx}=delay;
    Uscf_BP_n(idx,:)=interp1(vg,Uscf_mat(1,:,idx),vg_new)+min(min(E_all))-Ef;
end
vg_new=vg_new-v_off;
[r, theta]=meshgrid(vg_new,Angl);
XM=r.*cosd(theta);
YM=r.*sind(theta);
pcolor(XM,YM,curr_new/1000);
shading interp
hold on;
contour(XM,YM,curr_new/1000,1:8,'k','linewidth',1.5);
caxis([0 8]);
text(0.4,0.5,'n-type')
IV_data=curr_new([1 end],:)';
set(gca,'xtick',0:0.2:0.6,'ytick',0:0.2:0.6);
h=xlabel('V_G (V) [x-dir]');
set(h,'units','inches','position',[0.575 -0.2]);
h=ylabel('V_G (V) [y-dir]');
set(h,'units','inches','position',[-0.25 0.6]);
text(-0.12,0.53,'(a)');

axes('units','inches','position',[2.1 2.25 1.2 1.15]);
load data1/BP_p_angle_dense Vg Angl I NQ
NQ=squeeze(NQ);
Vel_old=I./NQ;
vd=0.6;
vg=Vg;
curr_old=I;
curr_new=zeros(length(Angl),21);
Vel_new_p=zeros(length(Angl),21);
delay_p=cell(length(Angl),1);
for idx=1:length(Angl)
    curr=curr_old(:,idx);
    v_off=interp1(log10(curr),vg,-1);
    v_on=v_off+vd;
    vg_new=linspace(v_off,v_on,21);
    curr_new(idx,:)=(10.^interp1(vg,log10(curr),vg_new));
    vel=Vel_old(:,idx);
    Vel_new_p(idx,:)=(interp1(vg,vel,vg_new));
    NQ_tmp=NQ(:,idx);
    delay=zeros(11,2);
    for idx_vg=1:11
        dQ=NQ_tmp(idx_vg+[0 30]);
        dI=curr(idx_vg+[0 30]);
        delay(idx_vg,:)=[(dI(2)/dI(1)) (dQ(2)-dQ(1))/(dI(2)-dI(1))];
    end
    delay_p{idx}=delay;
end
vg_new=vg_new-v_off;
[r, theta]=meshgrid(vg_new,Angl);
XM=r.*cosd(theta);
YM=r.*sind(theta);
pcolor(XM,YM,curr_new/1000);
shading interp
hold on;
contour(XM,YM,curr_new/1000,1:8,'k','linewidth',1.5);
caxis([0 8]);
text(0.4,0.5,'p-type')
IV_data=[IV_data curr_new([1 end],:)'];
set(gca,'xtick',0:0.2:0.6,'ytick',0:0.2:0.6);
h=xlabel('V_G (V) [x-dir]');
set(h,'units','inches','position',[0.6 -0.2]);
h=colorbar;
set(h,'unit','inches','position',[3.4 2.2 0.2 1.2]);
h=ylabel('Current (mA/\mum)');
set(h,'unit','inches','position',[1.7 0.575],'rotation',-90);
text(-0.12,0.53,'(b)');

axes('units','inches','position',[0.6 0.45 3 1.4]);
semilogy(vg_new,IV_data(:,1),style{1},'linewidth',1);
hold on;
semilogy(vg_new,IV_data(:,3),style{3},'linewidth',1);
semilogy(vg_new,IV_data(:,2),style{2},'linewidth',1);
semilogy(vg_new,IV_data(:,4),style{4},'linewidth',1);
set(gca,'box','off')
xlim([0 vd]);
ylim([1e-1 1e4]);
xlabel('V_G (V)')
ylabel('Current (\muA/\mum)')
axes('units','inches','position',[0.6 0.45 3 1.4]);
plot(vg_new,IV_data(:,1)/1000,style{1},'linewidth',1);
hold on;
plot(vg_new,IV_data(:,3)/1000,style{3},'linewidth',1);
plot(vg_new,IV_data(:,2)/1000,style{2},'linewidth',1);
plot(vg_new,IV_data(:,4)/1000,style{4},'linewidth',1);
set(gca,'color','none','yaxislocation','right','xticklabel',[],...
    'box','off','xaxislocation','top')
xlim([0 vd]);
ylim([0 10]);
h=ylabel('Current (mA/\mum)');
set(h,'rotation',-90,'units','inches','position',[3.3 0.7])
h=legend('n x-dir','p x-dir','n y-dir','p y-dir');
set(h,'units','inches','position',[0.58 1.27 0.9 0.55],'box','off');
text(-0.1,7500,'(c)');
set(gcf,'paperpositionmode','auto');
print -dtiff -r300 Fig2.tif

figure('units','inches','position',[5 1 4 2]);
load BP_mono_dense_EC
kx_plot=kx_plot*(2*pi/4.5694e-10);
ky_plot=ky_plot*(2*pi/3.3255e-10);
axes('units','inches','position',[0.475 0.45 1.5 1.5]);
Emin=min(min(E_all));
pcolor(kx_plot/1e9,ky_plot/1e9,E_all-Emin);
shading interp
caxis([0 1]);
hold on;
contour(kx_plot/1e9,ky_plot/1e9,E_all-Emin,0.05,'w','linewidth',1.5);
xlabel('k_x (1/nm)');
ylabel('k_y (1/nm)');
text(-10.5,8,'(a)')
load BP_mono_dense_EV
kx_plot=kx_plot*(2*pi/4.5694e-10);
ky_plot=ky_plot*(2*pi/3.3255e-10);
axes('units','inches','position',[2.45 0.45 1.5 1.5]);
Emin=min(min(E_all));
pcolor(kx_plot/1e9,ky_plot/1e9,E_all-Emin);
shading interp
caxis([0 1]);
hold on;
contour(kx_plot/1e9,ky_plot/1e9,E_all-Emin,0.05,'w','linewidth',1.5);
xlabel('k_x (1/nm)');
ylabel('k_y (1/nm)');
text(-10.5,8,'(b)')
colormap(flipud(colormap('cool')));
set(gcf,'paperpositionmode','auto');
print -dtiff -r300 Fig1.tif

figure('units','inches','position',[1 5 4 2]);
axes('units','inches','position',[0.4 0.45 1.5 1.5]);
plot(vg_new,Vel_new_n(1,:)/1e4,style{1}(1:end-1),'linewidth',1);
hold on;
plot(vg_new,Vel_new_p(1,:)/1e4,style{3}(1:end-1),'linewidth',1);
plot(vg_new,Vel_new_n(end,:)/1e4,style{2}(1:end-1),'linewidth',1);
plot(vg_new,Vel_new_p(end,:)/1e4,style{4}(1:end-1),'linewidth',1);
xlim(vg_new([1 end]));
ylim([0 35]);
h=xlabel('V_G (V)');
set(h,'units','inches','position',[0.75 -0.25]);
h=ylabel('<\nu> (10^4m/s)');
set(h,'units','inches','position',[-0.2 0.75]);
text(-0.12,34,'(a)');
text(0.02,16,'x-dir');
text(0.45,8,'y-dir');

axes('units','inches','position',[0.75 1.5 0.6 0.4],'fontsize',8);
plot(Angl,Vel_new_n(:,end)/1e4,'b-');
hold on;
plot(Angl,Vel_new_p(:,end)/1e4,'r--');
xlim(Angl([1 end]));
ylim([5 25]);
set(gca,'xtick',0:10:90,'xticklabel',{'X','','','','','','','','','Y'})
h=xlabel('Trans. dir.');
set(h,'units','inches','position',[0.3 -0.1]);
h=ylabel('<\nu_{sat}>');
set(h,'units','inches','position',[-0.18 0.2]);
text(2,14,{'V_G=','0.6V'},'fontsize',8);

axes('units','inches','position',[2.5 0.45 1.3 1.5]);
delay=delay_n{1};
semilogx(delay(:,1),delay(:,2)*20e-9/1e-12,'b-','linewidth',1);
hold on;
delay=delay_p{1};
semilogx(delay(:,1),delay(:,2)*20e-9/1e-12,'r--','linewidth',1);
delay=delay_n{end};
semilogx(delay(:,1),delay(:,2)*20e-9/1e-12,'b-','linewidth',1);
delay=delay_p{end};
semilogx(delay(:,1),delay(:,2)*20e-9/1e-12,'r--','linewidth',1);
ylim([0.06 0.34]);
%xlim([2 8]);
h=xlabel('I_{ON}/I_{OFF}');
set(h,'units','inches','position',[0.65 -0.25]);
h=ylabel('\tau (ps@L_{CH}=20nm)');
set(h,'units','inches','position',[-0.3 0.75]);
text(1e1,0.33,'(b)');
text(1e3,0.2,'y-dir');
text(5e4,0.1,'x-dir');
h=legend('n-type','p-type');
set(h,'units','inches','position',[2.9 1.3 0.7 0.4],'box','off');
set(gcf,'paperpositionmode','auto');
print -dtiff -r300 Fig3.tif

vtmp=linspace(0,0.6,31);
load ZP_data/Imonx.mat
itmp=interp1(vtmp,Imonx,vg_new);
IV_data_MoS2_3=itmp';
load ZP_data/Imony.mat
itmp=interp1(vtmp,Imony,vg_new);
IV_data_MoS2_3=[IV_data_MoS2_3 itmp'];
load ZP_data/Imopx.mat
itmp=interp1(vtmp,Imopx,vg_new);
IV_data_MoS2_3=[IV_data_MoS2_3 itmp'];
load ZP_data/Imopy.mat
itmp=interp1(vtmp,Imopy,vg_new);
IV_data_MoS2_3=[IV_data_MoS2_3 itmp'];

load ZP_data/Isinx.mat
itmp=interp1(vtmp,Isinx,vg_new);
IV_data_Si=itmp';
load ZP_data/Isiny.mat
itmp=interp1(vtmp,Isiny,vg_new);
IV_data_Si=[IV_data_Si itmp'];
load ZP_data/Isip.mat
itmp=interp1(vtmp,Isip,vg_new);
IV_data_Si=[IV_data_Si itmp'];

figure('units','inches','position',[5 5 4 2]);
axes('units','inches','position',[0.42 0.45 1.4 1.3]);
semilogy(vg_new,IV_data(:,1),style{1},'linewidth',1);
hold on;
semilogy(vg_new,IV_data_MoS2_3(:,1),style{2},'linewidth',1);
semilogy(vg_new,IV_data_Si(:,1),style{5},'linewidth',1);
semilogy(vg_new,IV_data(:,2),style{3},'linewidth',1);
semilogy(vg_new,IV_data_MoS2_3(:,2),style{4},'linewidth',1);
%semilogy(vg_new,IV_data_Si(:,2),style{6},'linewidth',1);
xlim(vg_new([1 end]));
ylim([1e-1 1e4]);
xlabel('V_G (V)')
h=ylabel('Current (\muA/\mum)');
set(h,'units','inches','position',[-0.25 0.65]);
set(gca,'box','off','ytick',[1e0 1e2 1e4])
text(0.02,5e3,'n-type');
text(-0.1,2e3,'(a)');
axes('units','inches','position',[0.42 0.45 1.4 1.3]);
plot(vg_new,IV_data(:,1)/1e3,style{1},'linewidth',1);
hold on;
plot(vg_new,IV_data_MoS2_3(:,1)/1e3,style{2},'linewidth',1);
plot(vg_new,IV_data_Si(:,1)/1e3,style{5},'linewidth',1);
plot(vg_new,IV_data(:,2)/1e3,style{3},'linewidth',1);
plot(vg_new,IV_data_MoS2_3(:,2)/1e3,style{4},'linewidth',1);
%plot(vg_new,IV_data_Si(:,2)/1e3,style{6},'linewidth',1);
xlim(vg_new([1 end]));
ylim([0 9]);
set(gca,'color','none','yaxislocation','right','xticklabel',[],...
    'box','off','xaxislocation','top')
%hl=legend('BP(x)','MoS_2(x)','Si','BP(y)','MoS_2(y)');
%set(hl,'units','inches','position',[1 1.75 2 0.2])%,...
    %'fontsize',8,'orientation','horizontal','box','on');

axes('units','inches','position',[2.3 0.45 1.4 1.3]);
semilogy(-vg_new,IV_data(:,3),style{1},'linewidth',1);
hold on;
semilogy(-vg_new,IV_data_MoS2_3(:,3),style{2},'linewidth',1);
semilogy(-vg_new,IV_data_Si(:,3),style{5},'linewidth',1);
semilogy(-vg_new,IV_data(:,4),style{3},'linewidth',1);
semilogy(-vg_new,IV_data_MoS2_3(:,4),style{4},'linewidth',1);
xlim(-vg_new([end 1]));
ylim([1e-1 1e4]);
xlabel('V_G (V)')
h=ylabel('Current (mA/\mum)');
set(h,'units','inches','position',[-0.3 0.65]);
set(gca,'box','off','ytick',[1e0 1e2 1e4],'yaxislocation','right',...
    'xtick',-0.6:0.2:0)
text(-0.02,5e3,'p-type','horizontalalignment','right');
text(-0.75,5e3,'(b)');

axes('units','inches','position',[2.3 0.45 1.4 1.3]);
plot(-vg_new,IV_data(:,3)/1e3,style{1},'linewidth',1);
hold on;
plot(-vg_new,IV_data_MoS2_3(:,3)/1e3,style{2},'linewidth',1);
plot(-vg_new,IV_data_Si(:,3)/1e3,style{5},'linewidth',1);
plot(-vg_new,IV_data(:,4)/1e3,style{3},'linewidth',1);
plot(-vg_new,IV_data_MoS2_3(:,4)/1e3,style{4},'linewidth',1);
xlim(-vg_new([end 1]));
ylim([0 9]);
set(gca,'color','none','xticklabel',[],...
    'box','off','xaxislocation','top')
set(gcf,'paperpositionmode','auto');
print -dtiff -r300 Fig4.tif