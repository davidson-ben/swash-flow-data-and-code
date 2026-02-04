%==========================================================================
% dye_image.m
%
% Visualize dye clouds.
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
%path = strcat('/Volumes/npujara/Ben Davidson/Queens 2024 Data/',Surface,'_A',ACT(1:2),'_C123/',ACT,'_',Surface,'/');
%load(strcat(path,'comb_xs_',ACT,'_',ref,'.mat'),'t_shore','x_shore')

[frame_time,time_cyc] = ensemble_avg_frames_dye;
    %frame_time is the time in seconds since the camera starts
    %time_cyc is the corresponding time in the swash cycle for the respective frame

% Starting Frames from images
start_frame = [848 1811 2601 3509 910 1878 2661];
end_frame =   [881 1846 2632 3541 932-2 1900 2675]; %ending before next wave

dye_drop = {'1','2','3','4','1','2','3','4'};

tests = 1:7;

for j = 1:length(tests)
    i = tests(j);
    %[xs, t_norm, bcl] = normalize_time_fn(ACT(1:2),'Smooth');
    load(strcat('../../7_dye/A00_Smooth',dye_drop{i},'.mat'))
    time = index;
    st = find(time==start_frame(i));
    ed = find(time==end_frame(i));
    time = time(st:ed);
    varX = varX(st:ed);
    t_n{i} = time_cyc(time);
    sdX{i} = sqrt(varX);
    x0(j) = x(st);
    t0(j) = t_n{i}(1);
end

%% Plot dye images


f2 = figure(2);
%f2.Position = [152    277   742   589];
clf

    I = imread("../../7_dye/dye_cloud_1.tif");
    
    imshow(I(:,:,1)) %plot red band
    clim([50 150]) %adjust clims for saturation

    %Add lines for dye spread
    hold on

    %plot var
    i = 1;
    load(strcat('../../7_dye/A00_Smooth',dye_drop{i},'.mat'))
    time = index;
    st = find(time==start_frame(i));
    ed = find(time==end_frame(i));
    time = time(st:ed);
    varX = varX(st:ed);
    varY = varY(st:ed);
    xo = x(st:ed);
    x = abs(x(st:ed)-7.7631)*(1000/1)*(1/mm_px);
    yo = y(st:ed);
    y = (y(st:ed)-9.7324)*(1000/1)*(1/mm_px);

    t_n = time_cyc(time);
    
    %alpha = fliplr(logspace(-0.5,0,length(x)));
    for j = [1 round(length(x)/2) length(x)]
        %plot(x(j),y(j),'ro','MarkerSize',10,'LineWidth',1)
        %alpha = 0.5;
        %c = [RGB(i,:) alpha(j)];
        
        sdx = sqrt(varX(j))*(1000/1)*(1/mm_px);
        plot([x(j)-2*sdx x(j)+2*sdx],[y(j) y(j)],'r-','linewidth',1.5)
        sdy = sqrt(varY(j))*(1000/1)*(1/mm_px);
        plot([x(j) x(j)],[y(j)-2*sdy y(j)+2*sdy],'r-','linewidth',1.5)
    end

    plot([xs0 xs0],[0 2056],'r--','linewidth',2)

    quiver(xs0,2056-175,-300,0,0,'r','linewidth',2,'MaxHeadSize',1)
    text(xs0-450,2056-75,'$x$','Color','red','Interpreter','latex','fontsize',50,'rotation',00)

    
