%==========================================================================
% figure_point_measurements.m
%
% Figure 2
% Plot point measurements and corresponding model for velocity and depth
% with alpha calculation.
% 
% Author: B. Davidson
% Last Updated: 3 October 2025
%==========================================================================

clear;
close all;
clc;


%% Load Data
load WG_locs.mat WG_locs %WG LOCATIONS

WG_locs_m = [-1.8 -1.4 -1 -0.7 -0.3 0.3];

load('A00_SWL.mat')
xs0 = X_swl;
clear X_swl;

load('cc_smooth.mat') %mm/px conversion


%% Add Antuono
Case = "03"; 
path_comb = strcat("R:\Ben Davidson\Queens 2024 Data\S_A00_C123\00",Case,"_comb_S");

if strcmp(Case,"01")
    alpha_measure = 2.4714 %2.3221 %2.3195;
    alpha_runup = 2.3221
    alp_m = 1.84;
    alp_r = 1.4666;
elseif strcmp(Case,"02")
    alpha_measure = 2.4886 %2.3401 %2.3436;
    alpha_runup = 2.3401
    alp_m = 1.88
    alp_r = 1.5146;
elseif strcmp(Case,"03")
    alpha_measure = 2.5419 %2.3765 %2.3936;
    alpha_runup = 2.3765;
    alp_m = 2.001;
    alp_r = 1.60889;
end


%load antuono solution - measured alpha
load(strcat('R:\Ben Davidson\0_Active_Processing\14 - Antuono\BD Code\alpha2_',string(alpha_measure),'.mat'),'X','T','D','xm','U','t0')
Xm{1} = X*1.8;
xmm{1} = xm*1.8;
Tm{1} = (T-t0)*1.3546;
Dm{1} = D*0.18;
Um{1} = U*1.3288;
clear X xm T D U t0

% %load antuono solution - run-up
% load(strcat('R:\Ben Davidson\0_Active_Processing\14 - Antuono\BD Code\alpha2_',string(alpha_runup),'.mat'),'X','T','D','xm','U','t0')
% Xm{2} = X*1.8;
% xmm{2} = xm*1.8;
% Tm{2} = (T-t0)*1.3546;
% Dm{2} = D*0.18;
% Um{2} = U*1.3288;
% clear X xm T D U t0




%% Sensors
x = [0  3.5;
     0    2.5;
     -0.5 2;
     -0.5 2;
     -1   1.5;
     -1   1.5;]; %x range of each plot

x = [-1 4; -1 4; -1 4; -1 4; -1 4; -1 4];

%xmodel = (xs0-WG_locs(:,1))*mm_px*(1/1000); %x location of each sensor
%xmodel(end) = 0.3;
xmodel = WG_locs_m;
%xmodel = flipud(xmodel);




%%
figure(1)
clf
tiledlayout(3,6,'TileSpacing','tight')
%% Depth
load(strcat(path_comb,"\combined_wave_gauges.mat"),'deA','deB','te')


%error('LOAD COMBINED DEPTH DATA')
load sw.mat %still water depth of sensors - maybe check on the other trials?
%load("../../Results/depth_ensemble_average.mat")
%return
% deB;   %ensemble average
% dptsB; %all points
% te_dB;    %ensemble average time
% time_cyc;          %all points time
edge = [1 0 0 0 0 0];
tl = {'F','E','D','C','B','A'};

%deB = flip(deB);

for i = 1:6
    nexttile
    %plot_depth_ref(dptsB{i},time_cyc,deB{i},te_dB{i},edge(i),tl{i},x(i,1),x(i,2),1,1,xmodel(i),xmm,Tm,Dm);
    %plot_depth_ref_2(deB{i},te,edge(i),tl{i},x(i,1),x(i,2),1,1,xmodel(i),xmm,Tm,Dm);
    plot_depth_ref_3(deB{i},te,edge(i),tl{i},x(i,1),x(i,2),1,1,xmodel(i),xmm,Tm,Dm)
    %plot_depth_model()
    %plot_depth_model()
    yline(0,'k-')
end


%% Velocity
load(strcat(path_comb,"\combined_vectrinos.mat"),'vel_out','te')

%% ADD ptv data
load(strcat(path_comb,'/ens_av_ptv_comb.mat'),'xp','tp','up')

%x = WG_locs(:,1); %x location of sensors in cross-shore (pixels)
%x_rel_swl = (xs0-x)*mm_px/1000;
%x_swl_round = round(x_rel_swl,2);

%x_sens = (xs0 - WG_locs([5 6],1))*mm_px/1000;
%x_sens(end) = 0.3;
x_sens = WG_locs_m(5:6);
%x_sens = x_rel_swl(5:6)

vE = nan(size(te));
vF = nan(size(te));
tw = diff(te(1:2))/2;
xw = 0.05/2;

for i = 1:length(te)
    tl = te(i) - tw;
    th = te(i) + tw;

    %valid time range
    if tl>=0 && th<=2
        t_val = tp>=tl & tp<=th;
    elseif tl<0
        t_val = tp>=(tl+2) | tp<=th;
    elseif th>2
        t_val = tp>=tl | tp<=(th-2);
    end

    %location E
    e_val = t_val;
    xh = x_sens(1) + xw;
    xl = x_sens(1) - xw;
    e_val(xp<xl) = 0;
    e_val(xp>xh) = 0;
    if sum(e_val)>=5
        vE(i) = median(up(e_val));
    else
        vE(i) = nan;
    end
    

    %location F
    f_val = t_val;
    xh = x_sens(2) + xw;
    xl = x_sens(2) - xw;
    f_val(xp<xl) = 0;
    f_val(xp>xh) = 0;
    if sum(f_val)>=5
        vF(i) = median(up(f_val));
    else
        vF(i) = nan;
    end
    
end


vel_out{5} = [vE vE vE];
vel_out{6} = [vF vF vF];

% figure(2)
% clf
% scatter(tp,xp,'k.')
% hold on
% x_sens = [-0.55 0.05];
% yline(x_sens)
% ylabel('x [m]')
% xlabel('t [s]')
% scatter(tp(f_val),xp(f_val),'r.')




%%



%nexttile
% V_B6_ea;   %ensemble average
% V_B6_trim; %all points
% t_V_ea;    %ensemble average time
% t;          %all points time

edge = [1 0 0 0 0 0];

for i = 1:6
    nexttile
    %plot_velocity_ref_2(vel_out{i},te,edge(i),x(i,1),x(i,2),1,1,xmodel(i),xmm,Tm,Um);
    plot_velocity_ref_3(vel_out{i},te,edge(i),x(i,1),x(i,2),1,1,xmodel(i),xmm,Tm,Um);
end




%% Alpha

alp{1} = alp_m*9.81*0.1; %9.81*0.1*2;
%alp{2} = alp_r*9.81*0.1; 

for i = 1:6
    nexttile
    plot_alpha_ref_3(vel_out{i},te,deB{i},edge(i),x(i,1),x(i,2),1,1,alp);
    tt = -1:0.1:4;
    gst = 9.81*0.1*tt;
    a_sw = gst + 2*sqrt(9.81*sw(i));
    plot(tt,gst,'b--','linewidth',2)
    %plot(tt,a_sw,'g--','linewidth',2)'

    %plot([-1 4],[alp alp],'g--')
    %plot([-1 4],[alp alp]+9.81*0.1*2,'g:')
end


exportgraphics(gcf,strcat("./save_figs/sensors_",Case,".pdf"),'ContentType','vector')