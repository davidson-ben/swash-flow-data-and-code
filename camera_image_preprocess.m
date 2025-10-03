%==========================================================================
% camera_image_preprocess.m

% Camera image pre-processing and calibration. 
% 
% Pre-process images from JAI camera images for SPTV, shoreline tracking,
% and dye dispersion experiments.
%
% Code removes lens distortion, converts to monochrome, inverts image, and
% saves the new image with a new name.
% 
% Author: B. Davidson
% Last Updated: 16 September 2025
%==========================================================================

clear;
close all;
clc;

% Camera calibration toolbox for MATLAB, J.-Y. Bouguet (2003)
addpath \Documents\MATLAB\TOOLBOX_calib  %download calibration toolbox and add the respective path here

% Original Camera Images
projectdir = strcat("Path to original camera images.");

% Save Folders
savedir = strcat("Path to save rectified images.");

maindir = pwd;

files = dir(strcat(projectdir,'/*.tiff')); %All images to be processed

cd(maindir)

load Calibration_Results_072424.mat fc alpha_c cc kc nx ny %load calibration results (Calibration results updated: 7/3/24)
KK = [fc(1) alpha_c*fc(1) cc(1);0 fc(2) cc(2) ; 0 0 1];
cd(projectdir)

% rename files
rename_start = str2double(string(regexp(files(1).name, '(\d{4})_','tokens')));
rename_end = str2double(string(regexp(files(end).name, '(\d{4})_','tokens')));

parfor id = 1:length(files) %loop through each file
    I = imread(files(id).name); %read image
    I_bw = double(mean(I,3)); %average to monochrome

    % UNDISTORT THE IMAGE:
    fprintf(1,strcat("Computing the undistorted image: ",string(id),"..."))
    [I_rect] = rect(I_bw(:,:,1),eye(3),fc,cc,kc,alpha_c,KK); %rectify via Bouguet (2003)
    I2 = ones(ny, nx,1);
    I2(:,:,1) = I_rect;
    fprintf(1,'done\n')
    newname = strcat(string(regexp(files(id).name, '(\d{4})_','tokens')),"_rect.tif"); %rename
    imwrite(uint8(I2), fullfile(savedir, newname)) %save image
end

t_start = 360; %camera starts 6 minutes after wave trigger
dt = 0.033; %time between camera frames from camera framerate (30.3 fps)
t_end = t_start + dt*(4000);
time = t_start:dt:t_end; %full time vecotor for camera
save(fullfile(savedir,'time_trim.mat'), 'time')
disp("DONE")