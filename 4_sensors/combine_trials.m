clear;
close all;
clc;

%% Setup
% All for wave angle = 0 degrees
Case = '02'; %Define wave case

path_t1 = strcat("R:\Ben Davidson\Queens 2024 Data\S_A00_C123\00",Case,"01_S");
path_t2 = strcat("R:\Ben Davidson\Queens 2024 Data\S_A00_C123\00",Case,"02_S");

%% PTV Field

%Check if necessary files exist for each trial
file = "ptv_xyuv.mat";
er(1) = check_file(path_t1,file);
er(2) = check_file(path_t2,file);
file = "time_trim.mat";
er(3) = check_file(path_t1,file);
er(4) = check_file(path_t2,file);
file = "shoreline.mat";
er(5) = check_file(path_t1,file);
er(6) = check_file(path_t2,file);

%confirm that cc_smooth and A00_SWL.mat both exist
file = "cc_smooth.mat";
path_cc = "R:\Ben Davidson\Queens 2024 Code\Image Preprocessing";
er(7) = check_file(path_cc,file);
file = "A00_SWL.mat";
path_swl = "R:\Ben Davidson\Queens 2024 Data\background";
er(8) = check_file(path_swl,file);

if sum(er)==0
    clear er
else
    error("SOME FILES MISSING")
end

%Load files TRIAL 01
load(strcat(path_t1,'\','ptv_xyuv.mat'))
load(strcat(path_t1,'\','time_trim.mat'))
load(strcat(path_t1,'\','shoreline.mat'))
%Load files both
load(strcat(path_cc,'\cc_smooth.mat'))
load(strcat(path_swl,'\A00_SWL.mat'))
x_ref = X_swl;



%Process Trial 01
[frame_time,time_cyc] = ensemble_avg_frames(x_ref, bins, wave);
    %frame_time is the time in seconds since the camera starts
    %time_cyc is the corresponding time in the swash cycle for the respective frame
[xp_t1, tp_t1, up_t1] = ea_uxt(time, frame_time, mm_px, wave, ptv_xyuv, time_cyc, x_ref);


%Load files TRIAL 02
load(strcat(path_t2,'\','ptv_xyuv.mat'))
load(strcat(path_t2,'\','time_trim.mat'))
load(strcat(path_t2,'\','shoreline.mat'))

%Process Trial 02
[frame_time,time_cyc] = ensemble_avg_frames(x_ref, bins, wave);
    %frame_time is the time in seconds since the camera starts
    %time_cyc is the corresponding time in the swash cycle for the respective frame
[xp_t2, tp_t2, up_t2] = ea_uxt(time, frame_time, mm_px, wave, ptv_xyuv, time_cyc, x_ref);

%Combine PTV field for trials 01 and 02
xp = [xp_t1; xp_t2];
tp = [tp_t1; tp_t2];
up = [up_t1; up_t2];

%%
figure(3)
plot(tp_t1,xp_t1,'r.')
hold on
plot(tp_t2,xp_t2,'k.')

save(strcat("R:\Ben Davidson\Queens 2024 Data\S_A00_C123\00",Case,"_comb_S\ens_av_ptv_comb.mat"),'xp','tp','up')