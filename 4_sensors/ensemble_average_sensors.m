%==========================================================================
% ensemble_average_sensors.m
%
% Ensemble average depth and velocity data from sensors and SPTV point
% measurements.
% 
% Author: B. Davidson
% Last Updated: 3 October 2025
%==========================================================================

clear;
close all;
clc;


path = "../Results";


% %Load data matrix
% %Load data matrix from original file: all packed up in one variable and hard to load   %load('/Volumes/npujara/Ben Davidson/0_Active_Processing/000 - QU Clean EA/ea_cycles.m')
% load('data_matrix_unpacked_Smooth_A00_C03.mat','VECTRINOB1') %unpacked file - can load individual variables
%     %sec_since_start
%     %time_local
%     %UDM(A/B)#
%     %Vectrino(A/B)#
%     %Wisc_SLB5
% 
% return


%load wave gauges
load('datalogger.mat','datalogger','dltime','column_desciption_datalogger')

%% ensemble average from shoreline and reference




[frame_time,time_cyc] = ensemble_avg_frames;
    %frame_time is the time in seconds since the camera starts
    %time_cyc is the corresponding time in the swash cycle for the respective frame

%% Trigger
%loading datalogger trigger and time
t_dl = dltime; %datetime string for each entry -- we will use this to sync velocimeter data
trig = datalogger(:,11); %trigger at each entry

% build time vector for for datalogger relative to camera start
    % use trigger to get start and end time if possible
    re = find(trig>4.5,1,'first'); %rising edge
    if re == 1
        re = nan;
    end
    fe = find(trig>4.5,1,'last'); %falling edge
    if fe == length(trig)
        fe = nan;
    end

    if isnan(re) && isnan(fe)
        error('rising edge and falling edge both nan')
    elseif isnan(re) && ~isnan(fe)
        disp('only use falling edge')
        t_end = t_dl(fe);
        t_cam_start = t_end - minutes(5); %time camera starts (5 minutes before trigger end)

        t_new = seconds(t_dl - t_cam_start); %time (in seconds) since camera start
    elseif ~isnan(re) && isnan(fe)
        disp('only use rising edge')
        error('still need to code')
    elseif ~isnan(re) && ~isnan(fe)
        disp('use both rising and falling edge')
        error('still need to code')
    end

%% Ensemble Average Depths

udmsB = {'B6','B5','B4','B3','B2','B1'};
udmsA = {'A6','A5','A3'};

udmsB_inds = [7 6 5 4 3 2];
udmsA_inds = [10 9 8];

for i = 1:length(udmsB)
    depth_dl = datalogger(:,udmsB_inds(i)); %depth at each entry
    [deB{i},te_dB{i},dptsB{i}] = ens_av(depth_dl,t_new,frame_time,time_cyc,5); %de is depth [Q1 Q2 Q3] at te
end

for i = 1:length(udmsA)
    depth_dl = datalogger(:,udmsA_inds(i)); %depth at each entry
    [deA{i},te_dA{i},dptsA{i}] = ens_av(depth_dl,t_new,frame_time,time_cyc,5); %de is depth [Q1 Q2 Q3] at te
end


save(strcat(path,'/depth_ensemble_average.mat'),'deA','te_dA','dptsA','deB','te_dB','dptsB','time_cyc')

%% Ensemble Average Velocities

%save and just load all velocimeter data
%loop through to solve all sensors

%same time for all vectrinos
load('Velocimeters_Smooth_A00_C03_T01.mat','time_local')
tu = seconds(time_local - t_cam_start);

%B6
ueB{1} = [];
te_uB{1} = [];
uptsB{1} = [];

%B5
ueB{2} = [];
te_uB{2} = [];
uptsB{2} = [];

%B4
load('Velocimeters_Smooth_A00_C03_T01.mat','VECTRINOB4')
u = VECTRINOB4.gVelX(end,:);
[ueB{3},te_uB{3},uptsB{3}] = ens_av(u,tu,frame_time,time_cyc,5);

%B3
load('Velocimeters_Smooth_A00_C03_T01.mat','VECTRINOB3')
u = VECTRINOB3.gVelX(end,:);
[ueB{4},te_uB{4},uptsB{4}] = ens_av(u,tu,frame_time,time_cyc,5);

%B2
load('Velocimeters_Smooth_A00_C03_T01.mat','VECTRINOB2')
u = VECTRINOB2.gVelX(end,:);
[ueB{5},te_uB{5},uptsB{5}] = ens_av(u,tu,frame_time,time_cyc,5);

%B1
load('Velocimeters_Smooth_A00_C03_T01.mat','VECTRINOB1')
u = VECTRINOB1.gVelX(end,:);
[ueB{6},te_uB{6},uptsB{6}] = ens_av(u,tu,frame_time,time_cyc,5);


%% Ensemble Avgerage PTV

%load PTV data (by frame)

%load track details
load("../Results/track_details_px_f.mat") %all particles position and data
load("../Results/time_trim.mat") %time
%B is each particle

N = length(time)-1;
x = nan(N,length(B));
y = x;
u = x;
v = x;

%loop through each particle and sort data
for p = 1:length(B)
    part = B{p};
    for i = 1:size(part,1)
        if all(~isnan(part(i,5:6))) && part(i,11) == 0 %if velocities are not nan and particle is not beached
            %ptv_xyuv{part(i,3)}(end+1,:) = [part(i,[1 2 5 6])];%append particle location and velocity to the ptv frame

            t = part(i,3);

            x(t,p) = part(i,1);
            y(t,p) = part(i,2);
            u(t,p) = part(i,5);
            v(t,p) = part(i,6);
        end
    end
end

%% Individual Points
load WG_locs.mat B6_px B5_px
%B6
x_b6 = B6_px(1);
y_b6 = B6_px(2);

validX = nan(size(x));

w = 150; %half width [px]

validX(x<x_b6+w & x>x_b6-w) = 1;

uvalid = u.*validX;

u_b6 = mean(uvalid,2,'omitnan');

val = sum(validX,2,'omitnan');

thresh = 1;

u_b6(val<thresh | isnan(val)) = nan; %cross shore velocity [px/frame]
load('../1_camera_preprocessing/cc_smooth.mat')


u_b6 = -1*u_b6*30.3*mm_px*(1/1000);
t_b6 = time(1:end-1);

[ueB{1},te_uB{1},uptsB{1}] = ens_av(u_b6,t_b6-360,frame_time,time_cyc,10);

%B5
x_b5 = B5_px(1);
y_b5 = B5_px(2);

validX = nan(size(x));

w = 150; %half width [px]

validX(x<x_b5+w & x>x_b5-w) = 1;

uvalid = u.*validX;

u_b5 = mean(uvalid,2,'omitnan');

val = sum(validX,2,'omitnan');

thresh = 1;

u_b5(val<thresh | isnan(val)) = nan; %cross shore velocity [px/frame]


u_b5 = -1*u_b5*30.3*mm_px*(1/1000);
t_b5 = time(1:end-1);

[ueB{2},te_uB{2},uptsB{2}] = ens_av(u_b5,t_b5-360,frame_time,time_cyc,10);


%save ensemble averaged velocity data
save(strcat(path,'/velocity_ensemble_average.mat'),'ueB','te_uB','uptsB')









function [xe,te,x_ft] = ens_av(x,t,frame_time,time_cyc,n_thresh)
    %x is the value to be ensemble averaged measured on time t
    %t is the time (s) for the value to be ensemble averaged, where 0 is the camera start
    
    %interpolate to frame time
    x_ft = interp1(t,x,frame_time);

    te = linspace(0,2,61)'; %ensemble average time (0 and 2 are the same)
    xe = nan(length(te),3); %ensemble averaged measured value (empty array) [Q1 Q2 Q3]

    dt = diff(te(1:2)); %te time-steps

    for i = 1:size(xe,1)
        %window of dt - 1/2 on either side of the time indexed
        t_low = te(i)-dt/2;
        t_high = te(i)+dt/2;

        %since swash is cycle, we need to wrap around on the ends
        %define indices where time is within the window (inds)
        if t_low < 0
            t_low = t_low + 2;
            inds = find(time_cyc < t_high | time_cyc > t_low);
        elseif t_high > 2
            t_high = t_high - 2;
            inds = find(time_cyc < t_high | time_cyc > t_low);
        else
            inds = find(time_cyc < t_high & time_cyc > t_low);
        end
    
        not_nan = sum(~isnan(x_ft(inds))); %number of non-nan entries
        if not_nan >= n_thresh %at least n_thresh non-nan values within a window requred to solve for ensemble value, otherwise left as nan
            xe(i,2) = median(x_ft(inds),'omitnan');
            xe(i,1) = quantile(x_ft(inds),0.25);
            xe(i,3) = quantile(x_ft(inds),0.75);
        end
    end

end