%==========================================================================
% dye_cloud_4bcde.m
%
% Figure 4b-e
% Dye cloud in x-t space with shoreline trajectory.  Dye cloud is shaded
% corresponding to specific dispersion regimes (swash stretching and
% enhanced).
% 
% Author: B. Davidson
% Last Updated: 3 October 2025
%==========================================================================

clear;
close all;
clc;

load("../../1_camera_preprocessing/A00_SWL.mat")
xs0 = X_swl;
clear X_swl;

load('../../1_camera_preprocessing/cc_smooth.mat') %mm/px conversion
path = "../../7_dye";

load("../../Results/combined_shoreline_ea.mat",'t_shore','x_shore')


[frame_time,time_cyc] = ensemble_avg_frames_dye;
    %frame_time is the time in seconds since the camera starts
    %time_cyc is the corresponding time in the swash cycle for the respective frame


% Starting Frames from images
start_frame = [848 1811 2601 3509];%[848 1811 2601 3509];
start_2nd = [910 1878 2661 nan];
end_frame =   [930 1900 2675 3548]; %ending before next wave

dye_drop = {'1','2','3','4','1','2','3','4'};

close all;
figure(3)
clf
tiledlayout(2,2,'TileSpacing','tight')

%% BLUE

%plot based on start of each cycle
[t, t_2nd, xc, sdX] = load_cloud(1,dye_drop,start_frame,start_2nd,end_frame,time_cyc,mm_px,xs0);
setup_plot(0,1);
plot_cloud(t,xc,sdX)

t_crit = [1.62 2.93 3.28 1.97];

%Swash stretching 1
t0 = t(1);
t1 = t_crit(1)*t0;
plot_swashstretching(t0,t1,t,xc,sdX)

%Swash stretching 2
t2 = t_crit(2)*t0;
t3 = t_crit(3)*t0;
plot_swashstretching(t2,t3,t,xc,sdX)

%Swash stretching 3
t0_2nd = t_2nd(1);
t4 = t_crit(4)*(t0_2nd-2)+2;
plot_swashstretching(t0_2nd,t4,t,xc,sdX)

%ENHANCED Stretching 1
te1 = t1;
te2 = t2;
plot_enhanced(t1,t2,t,xc,sdX)

%Plot shoreline (Combined from 000301)
plot(t_shore,x_shore,'k-','linewidth',2)
plot(t_shore+2,x_shore,'k-','linewidth',2)

%% Orange
%plot based on start of each cycle
[t, t_2nd, xc, sdX] = load_cloud(2,dye_drop,start_frame,start_2nd,end_frame,time_cyc,mm_px,xs0);
setup_plot(0,0);
plot_cloud(t,xc,sdX)

t_crit = [1.65 3.31 4.23];

%Swash stretching 1
t0 = t(1);
t1 = t_crit(1)*t0;
plot_swashstretching(t0,t1,t,xc,sdX)

%Swash stretching 2
t2 = t_crit(2)*t0;
t3 = t_crit(3)*t0;
plot_swashstretching(t2,t3,t,xc,sdX)

%ENHANCED Stretching 1
te1 = t1;
te2 = t2;
plot_enhanced(t1,t2,t,xc,sdX)

%Plot shoreline (Combined from 000301)
plot(t_shore,x_shore,'k-','linewidth',2)
plot(t_shore+2,x_shore,'k-','linewidth',2)

%% Yellow
%plot based on start of each cycle
[t, t_2nd, xc, sdX] = load_cloud(3,dye_drop,start_frame,start_2nd,end_frame,time_cyc,mm_px,xs0);
setup_plot(1,1);
plot_cloud(t,xc,sdX)

t_crit = [1.8 2.96 3.25];

%Swash stretching 1
t0 = t(1);
t1 = t_crit(1)*t0;
plot_swashstretching(t0,t1,t,xc,sdX)

%Swash stretching 2
t2 = t_crit(2)*t0;
t3 = t_crit(3)*t0;
plot_swashstretching(t2,t3,t,xc,sdX)

%ENHANCED Stretching 1
te1 = t1;
te2 = t2;
plot_enhanced(t1,t2,t,xc,sdX)

%Plot shoreline (Combined from 000301)
plot(t_shore,x_shore,'k-','linewidth',2)
plot(t_shore+2,x_shore,'k-','linewidth',2)


%% Purple
%plot based on start of each cycle
[t, t_2nd, xc, sdX] = load_cloud(4,dye_drop,start_frame,start_2nd,end_frame,time_cyc,mm_px,xs0);
setup_plot(1,0);
plot_cloud(t,xc,sdX)

t_crit = [1.65 2.54 2.90];

%Swash stretching 1
t0 = t(1);
t1 = t_crit(1)*t0;
plot_swashstretching(t0,t1,t,xc,sdX)

%Swash stretching 2
t2 = t_crit(2)*t0;
t3 = t_crit(3)*t0;
plot_swashstretching(t2,t3,t,xc,sdX)

%ENHANCED Stretching 1
te1 = t1;
te2 = t2;
plot_enhanced(t1,t2,t,xc,sdX)

%Plot shoreline (Combined from 000301)
plot(t_shore,x_shore,'k-','linewidth',2)
plot(t_shore+2,x_shore,'k-','linewidth',2)








function plot_cloud(t,xc,sdX)
    tp = [t fliplr(t)];
    xp = [xc-2*sdX fliplr(xc+2*sdX)];
    fill(tp,xp,'k','facealpha',0.3,'edgecolor','k')
end


function plot_swashstretching(t0,t1,t,xc,sdX)
ti = t(t<t1 & t>=t0);
xi = xc(t<t1 & t>=t0);
si = sdX(t<t1 & t>=t0);

    tp = [ti fliplr(ti)];
    xp = [xi-2*si fliplr(xi+2*si)];
    fill(tp,xp,'red','facealpha',0.7)
end

function plot_enhanced(t1,t2,t,xc,sdX)
    inds = find(t<t2 & t>t1);
    inds = [inds(1)-1 inds inds(end)+1];

    ti = t(inds);
    xi = xc(inds);
    si = sdX(inds);

    tp = [ti fliplr(ti)];
    xp = [xi-2*si fliplr(xi+2*si)];

    fill(tp,xp,'blue','facealpha',0.7)
end

function setup_plot(label_x, label_y)
    nexttile
    hold on
    set(gca,'FontSize',25)
    yticks(-0.5:0.5:1.0)
    ylim([-0.6 1.4])
    xticks(0:0.5:4)
    if label_x == 1
        xticklabels({'0','','1','','2','','3',''})
        xtickangle(0)
        xlabel('$t$ [s]','interpreter','latex')
    else
        xticklabels({})
    end

    if label_y == 1
        yticklabels({'-0.5', '0.0', '0.5', '1.0'})
        ylabel('$x$ [m]','interpreter','latex')
    else
        yticklabels({})
    end
    xlim([0 3.5])
    clim([-0.75 0.75])
    set(gca,'Tickdir','both')
    set(gca,'TickLabelInterpreter','latex')
end

function [t, t_2nd, xc, sdX] = load_cloud(i,dye_drop,start_frame,start_2nd,end_frame,time_cyc,mm_px,xs0)
    load(strcat('../../7_dye/A00_SMOOTH',dye_drop{i},'.mat'))
    time = index;
    st = find(time==start_frame(i));
    s2 = find(time==start_2nd(i)); %start 2nd cycle
    ed = find(time==end_frame(i));
    time_dye = time(st:ed);
    time_2nd = time(s2:ed);
    varX = varX(st:ed);
    varY = varY(st:ed);
    t = time_cyc(time_dye);
    t_2nd = time_cyc(time_2nd)+2;
    sdX = sqrt(varX);
    sdY = sqrt(varY);
    %xc = x(st:ed);
    xcp = abs(x(st:ed)-7.7631)*(1000/1)*(1/mm_px);
    xc = (xs0 - xcp)*mm_px*(1/1000);
    %yc = y(st:ed);
    %ycp = (y(st:ed)-9.7324)*(1000/1)*(1/mm_px);
    if i == 1 || i == 2 || i == 3
        t(find(t==min(t)):end) = t(find(t==min(t)):end)+2;
    end
end