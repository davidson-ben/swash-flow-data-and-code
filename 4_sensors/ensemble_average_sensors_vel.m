%==========================================================================
%==========================================================================

clear;
close all;
clc;

%% Setup
% All for wave angle = 0 degrees
Case = '03'; %Define wave case

%% Maybe add file check?

%% Load files
%Load files both
load('1_camera_preprocessing\cc_smooth.mat')
load('1_camera_preprocessing\A00_SWL.mat')
load(strcat("Results\Case","Case","\final_shoreline_ea_comb.mat"))

%Load files TRIAL 01
load(strcat("Results\Case",Case,'\','shoreline_t1.mat'))
wave1 = wave;
clear wave;
x_ref = X_swl;
%Process Trial 01
[frame_time1,time_cyc1] = ensemble_avg_frames(x_ref, bins, wave1);
    %frame_time is the time in seconds since the camera starts
    %time_cyc is the corresponding time in the swash cycle for the respective frame

%Load files TRIAL 02
load(strcat("Results\Case",Case,'\','shoreline_t2.mat'))
wave2 = wave;
clear wave;
%Process Trial 02
[frame_time2,time_cyc2] = ensemble_avg_frames(x_ref, bins, wave2);
    %frame_time is the time in seconds since the camera starts
    %time_cyc is the corresponding time in the swash cycle for the respective frame



%Load Sensors
%load(path_sensors,'VECTRINOB1','VECTRINOB2','VECTRINOB3','VECTRINOB4','WavePaddle_Voltage','UDMB1')
load(strcat("RAW_SMOOTH_Angle00_Case",Case,"R-2.mat"),'VECTRINOB1','VECTRINOB2','VECTRINOB3','VECTRINOB4','WavePaddle_Voltage','UDMB1')


%% Interpolate velocity data to datalogger time
%load data from trial 1
dltime1 = UDMB1.time_local{1};
trig_t1 = WavePaddle_Voltage{1};

V1234_t1(:,1) = interp1(VECTRINOB1.time_local{1},VECTRINOB1.VelX{1}(end,:),dltime1);
V1234_t1(:,2) = interp1(VECTRINOB2.time_local{1},VECTRINOB2.VelX{1}(end,:),dltime1);
V1234_t1(:,3) = interp1(VECTRINOB3.time_local{1},VECTRINOB3.VelX{1}(end,:),dltime1);
V1234_t1(:,4) = interp1(VECTRINOB4.time_local{1},VECTRINOB4.VelX{1}(end,:),dltime1);


%load data from trial 2
dltime2 = UDMB1.time_local{2};
trig_t2 = WavePaddle_Voltage{2};

V1234_t2(:,1) = interp1(VECTRINOB1.time_local{2},VECTRINOB1.VelX{2}(end,:),dltime2);
V1234_t2(:,2) = interp1(VECTRINOB2.time_local{2},VECTRINOB2.VelX{2}(end,:),dltime2);
V1234_t2(:,3) = interp1(VECTRINOB3.time_local{2},VECTRINOB3.VelX{2}(end,:),dltime2);
V1234_t2(:,4) = interp1(VECTRINOB4.time_local{2},VECTRINOB4.VelX{2}(end,:),dltime2);

%% Process

%trigger
t_dl{1} = dltime1;
t_dl{2} = dltime2;
trig{1} = trig_t1;
trig{2} = trig_t2;


vel{1} =V1234_t1;
vel{2} = V1234_t2;

frame_time{1} = frame_time1;
frame_time{2} = frame_time2;
time_cyc{1} = time_cyc1;
time_cyc{2} = time_cyc2;



%% Process Velocity
[vel_out, te] = ens_vectrino(t_dl,trig, vel, frame_time, time_cyc);



%%
pl = 1;
if pl == 1
    figure(1)
    clf
    hold on
    sB = 4;
    plot(te,vel_out{sB}(:,2),'r-','linewidth',1)
    plot(te,vel_out{sB}(:,[1 3]),'r--','linewidth',0.5)
end




function [vel_out, te] = ens_vectrino(t_dl,trig, vel, frame_time, time_cyc)
    
    %loading datalogger trigger and time
    %t_dl = dltime; %datetime string for each entry -- we will use this to sync velocimeter data
    %trig = datalogger(:,11); %trigger at each entry
    
    % build time vector for for datalogger relative to camera start
    % use trigger to get start and end time if possible
    for trial = 1:2
        re = find(trig{trial}>4.5,1,'first'); %rising edge
        if re == 1
            re = nan;
        end
        fe = find(trig{trial}>4.5,1,'last'); %falling edge
        if fe == length(trig{trial})
            fe = nan;
        end
    
        if isnan(re) && isnan(fe)
            %No rising nor falling edge - no trigger to compute time
            %relative to camera.
            error('rising edge and falling edge both nan')
        elseif isnan(re) && ~isnan(fe)
            %Only the falling edge is available and/or resolved.  Only use
            %falling edge to determine time relative to camera
            disp('only use falling edge')
            t_end = t_dl{trial}(fe);
            t_cam_start = t_end - minutes(5); %time camera starts (5 minutes before trigger end)
    
            t_new = seconds(t_dl{trial} - t_cam_start); %time (in seconds) since camera start
        elseif ~isnan(re) && isnan(fe)
            %Only the rising edge is available and/or resolved.  Only use
            %the rising edge to determine the time relative to the camera.
            disp('only use rising edge')

            t_start = t_dl{trial}(re);
            t_cam_start = t_start + minutes(6); %time camera starts (6 minutes after trigger start)
    
            t_new = seconds(t_dl{trial} - t_cam_start); %time (in seconds) since camera start
        elseif ~isnan(re) && ~isnan(fe)
            %Both rising and falling edge available and/or resolved.  Use
            %both rising and falling edge to determine the time relative to
            %the camera.
            disp('use both rising and falling edge')

            t_start = t_dl{trial}(re);
            t_end = t_dl{trial}(fe);
            t_cam_start_sttrig = t_start + minutes(6); %time camera starts (5 minutes before trigger end)
            t_cam_start_edtrig = t_end - minutes(5); %time camera starts (5 minutes before trigger end)

            t_cam_start = mean([t_cam_start_sttrig t_cam_start_edtrig]);
    
            t_new = seconds(t_dl{trial} - t_cam_start); %time (in seconds) since camera start
        end
        time{trial} = t_new;
    end

    %Results in time{trial} that uses all available triggers to best line
    %up the data.
    
    %% Ensemble Average Depths

    % A3 - ensemble average
    time; %time on the sensor signal relative to camera start
    frame_time; %cumulative time on the camera relative to camera start
    time_cyc; %cycle time on camera

    %interpolate sensor data to frame_time
    for i = 1:4 %across 4 vectrinos
        [vel_out{i}, te] = ens_av_two_trials(time, vel, i, frame_time, time_cyc, 5);
    end



    
    
    % for i = 1:length(udmsB)
    %     depth_dl = datalogger(:,udmsB_inds(i)); %depth at each entry
    %     [deB{i},te_dB{i},dptsB{i}] = ens_av(depth_dl,t_new,frame_time,time_cyc,5); %de is depth [Q1 Q2 Q3] at te
    % end
    % 
    % for i = 1:length(udmsA)
    %     depth_dl = datalogger(:,udmsA_inds(i)); %depth at each entry
    %     [deA{i},te_dA{i},dptsA{i}] = ens_av(depth_dl,t_new,frame_time,time_cyc,5); %de is depth [Q1 Q2 Q3] at te
    % end
end




function [xe,te] = ens_av_two_trials(t,x,sens,frame_time,time_cyc,n_thresh)
    %x is the value to be ensemble averaged measured on time t
    %t is the time (s) for the value to be ensemble averaged, where 0 is the camera start
    
    %interpolate to frame time
    x_ft{1} = interp1(t{1},x{1}(:,sens),frame_time{1});
    x_ft{2} = interp1(t{2},x{2}(:,sens),frame_time{2});

    te = linspace(0,2,61)'; %ensemble average time (0 and 2 are the same)
    xe = nan(length(te),3); %ensemble averaged measured value (empty array) [Q1 Q2 Q3]

    dt = diff(te(1:2)); %te time-steps

    for i = 1:size(xe,1)
        for trial = 1:2
            %window of dt - 1/2 on either side of the time indexed
            t_low = te(i)-dt/2;
            t_high = te(i)+dt/2;
    
            %since swash is cycle, we need to wrap around on the ends
            %define indices where time is within the window (inds)

            if t_low < 0
                t_low = t_low + 2;
                inds{trial} = find(time_cyc{trial} < t_high | time_cyc{trial} > t_low);
            elseif t_high > 2
                t_high = t_high - 2;
                inds{trial} = find(time_cyc{trial} < t_high | time_cyc{trial} > t_low);
            else
                inds{trial} = find(time_cyc{trial} < t_high & time_cyc{trial} > t_low);
            end
        end
    
        not_nan = sum(~isnan(x_ft{1}(inds{1}))) + sum(~isnan(x_ft{2}(inds{2}))); %number of non-nan entries
        if not_nan >= n_thresh %at least n_thresh non-nan values within a window requred to solve for ensemble value, otherwise left as nan
            xe(i,2) = median([x_ft{1}(inds{1}) x_ft{2}(inds{2})],'omitnan');
            xe(i,1) = quantile([x_ft{1}(inds{1}) x_ft{2}(inds{2})],0.25);
            xe(i,3) = quantile([x_ft{1}(inds{1}) x_ft{2}(inds{2})],0.75);
        end
    end

end


