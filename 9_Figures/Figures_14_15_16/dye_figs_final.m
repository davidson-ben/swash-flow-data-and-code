%==========================================================================
% dye_figs_final.m
%
% Figure 4a
% Plot dye dispersion of all clouds.
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

%%

RGB = orderedcolors("gem");
% 
% figure(2)
% clf
% hold on
% for i = 1:4
%     plot(t0(i),x0(i),'o','linewidth',2,'markersize',10,'color',RGB(i,:),'markerfacecolor',RGB(i,:))
% end
% title('Initial Conditions')
% xlabel('Cycle Time [s])')
% ylabel('Cloud Centroid [m]')
% set(gca,'FontSize',15)


figure(1)
clf
hold on

pf_d1 = 1.65;
x = linspace(2,4.5,10);
loglog(x,pf_d1*x.^(4/3),'k:','linewidth',2)

pf_d2 = 1.425;
loglog(x,pf_d2*x.^(4/3),'k:','linewidth',2)

pf_d3 = 1.27;
loglog(x,pf_d3*x.^(4/3),'k:','linewidth',2)

%pf_d4 = 1.7;
%loglog(x,pf_d4*x.^(4/3),'k:','linewidth',2)



for it = 1:4
    i = tests(it);
    t3 = (t_n{i})/t_n{i}(1);
    y3 = (sdX{i}).^2/sdX{i}(1)^2;

    loglog(t3,y3,'o-','Color',RGB(i,:),'MarkerFaceColor',RGB(i,:),'linewidth',1,'markersize',5)

end

for it = 5:7
    i = tests(it);
    t3 = (t_n{i})/t_n{i}(1);
    y3 = (sdX{i}).^2/sdX{i}(1)^2;

    loglog(t3,y3,'d-','Color',RGB(i-4,:),'linewidth',1,'markersize',5)

end




set(gca,'XScale','log','YScale','log')
set(gca,'FontSize',25)
xlabel('$t/t_0$','interpreter','latex')
%ylabel('$\sigma_x^2/\sigma_{x0}^2$','interpreter','latex')
ylabel('$ \frac{\sigma_x^2}{ \sigma_{x0}^2} $','interpreter','latex','rotation',0,'FontSize',30)
ylim([1 11])
yticks([1 2:2:10])
xticks([1 1.5 2 2.5 3 3.5 4])
xticklabels({'1.0', '1.5', '2.0', '2.5', '3.0', '3.5', '4.0'})
xlim([1 4.5])
set(gca,'Tickdir','both')
set(gca,'TickLabelInterpreter','latex')

%Swash Stretching
x = linspace(1,3,10);
pf_start = 1;
loglog(x,pf_start*x.^(4/3),'k--','linewidth',2)
txt = '$ \frac{\sigma_{x}^2}{\sigma_{x0}^2}  = \left( \frac{t}{t_0} \right) ^{4/3}$';
text(1.1,1.7,txt,'FontSize',25,'interpreter','latex','Rotation',30)




txt = '$\sim \left( \frac{t}{t_0} \right) ^{4/3}$';
text(2.0,6.0,txt,'FontSize',25,'interpreter','latex','rotation',30)


ylim([1 11])
yticks([1 2:2:10])
yticklabels([1 2 4 6 8])
%set(gca,"TickLabelInterpreter",'latex')
