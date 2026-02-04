%==========================================================================
% contours_data.m
%
% Figure 3 - data
% Plot individual panels for Figure 3 - data results.
% Run each section individually to create each panel.
% 
% Author: B. Davidson
% Last Updated: 3 October 2025
%==========================================================================

clear;
close all;
clc;

%%
%%%%%%%%%%%%
%DATA
%%%%%%%%%%%%

Case = "03";
path_comb = strcat("R:\Ben Davidson\Queens 2024 Data\S_A00_C123\00",Case,"_comb_S");

%load(strcat(path_comb,'/ea_depth_field.mat'),'Xd','Td','Dd')
load(strcat(path_comb,'/ens_av_ptv_comb.mat'))
load(strcat(path_comb,'/final_shoreline_ea_comb.mat'),'t_shore','x_shore')
load(strcat(path_comb,'/ea_depth_field_comb.mat'))



%% Solve Contours

%x windows
x_width = 0.1; %0.1 m
x_step = 0.05; %0.05 m
%start and end in the swash coords
x_start = -0.75;
x_end = 1.5;
%center of each x window
x_center = x_start:x_step:x_end;

%t windows
t_width = 0.2;
t_step = 0.1;
%start and end in swash time
t_start = 0;
t_end = 2;
%center of each t window
t_center = t_start:t_step:t_end;

[T,X] = meshgrid(t_center,x_center);
N = nan(length(x_center),length(t_center));
nthresh = 10;
U = nan(length(x_center),length(t_center));
D = nan(length(x_center),length(t_center));

for ti = 1:length(t_center)
    t_low = t_center(ti) - t_width/2;
    t_high = t_center(ti) + t_width/2;

    if t_low < 0
        validt = tp>(t_low+2) | tp<t_high;
    elseif t_high >2
        validt = tp>t_low | tp<(t_high-2);
    else
        validt = tp>t_low & tp<t_high;
    end
    for xi = 1:length(x_center)
        x_low = x_center(xi) - x_width/2;
        x_high = x_center(xi) + x_width/2;
        validx = xp>x_low & xp<x_high;

        uvalid = (validx==1 & validt==1);
        N(xi,ti) = sum(uvalid);
        if N(xi,ti)>nthresh
            U(xi,ti) = median(up(uvalid));
            
            D(xi,ti) = thiel_sen(xp(uvalid),up(uvalid));
            %s = polyfit(xp(uvalid),up(uvalid),1);
            %D(xi,ti) = s(1);
        end
    end
end

addpath('./crameri_colormaps')

if strcmp(Case,"01")
    for i = find(T(1,:)>0.3 & T(1,:)<0.7)
        tm = T(1,i)+2;
        shore_i = find(abs(t_shore-tm)==min(abs(t_shore-tm)));
        shore = x_shore(shore_i);
        if isnan(shore)
            shore = shore_keep;
        else
            shore_keep = shore;
        end
        rm = X(:,1)>shore;%entried ahead of shoreline to remove
        U(rm,i) = nan;
        D(rm,i) = nan;
    end
    for i = find(T(1,:)>0.7 & T(1,:)<1.25)
        tm = T(1,i);
        shore_i = find(abs(t_shore-tm)==min(abs(t_shore-tm)));
        shore = x_shore(shore_i);
        rm = X(:,1)>shore;%entried ahead of shoreline to remove
        U(rm,i) = nan;
        D(rm,i) = nan;
    end
    %U(,T(1,:)>0.7 & T(1,:)<1.25) = nan;
end

%% Solve Alpha

Dep_int = interp2(Td,Xd,Dd,T,X);

A = U + 2*sqrt(Dep_int*9.81) + 9.81*(1/10)*T;

%
% figure(3)
% clf
% plot(t_shore,x_shore,'-')
p = polyfit(t_shore(~isnan(x_shore)),x_shore(~isnan(x_shore)),2);
% hold on
% xp = -0.5:0.1:2.5;
% plot(xp,polyval(p,xp),'-')

%clean up alpha cylces
% for ti = 1:size(A,2)
%     if T(1,ti)>1
%         continue
%     end
% 
%     %get shoreline at this time
%     xs_q = polyval(p,T(1,ti));
%     A(X(:,ti)>xs_q,ti) = A(X(:,ti)>xs_q,ti) + 2;
% end

save(strcat(path_comb,'/plot_uhad.mat'))



%% Velocity Data Plot





figure(1)
clf
hold on
contourf(T-2,X,U,-0.75:0.1:0.75,'linestyle','none')
contourf(T,X,U,-0.75:0.1:0.75,'linecolor','none')
contourf(T+2,X,U,-0.75:0.1:0.75,'linecolor','none')
plot(t_shore-2,x_shore,'k.-','linewidth',2,'markersize',10)
plot(t_shore,x_shore,'k.-','linewidth',2,'markersize',10)
plot(t_shore+2,x_shore,'k.-','linewidth',2,'markersize',10)
set(gca,'FontSize',25)

%labels and ticks
xlabel('$t$ [s]','interpreter','latex','FontSize',40)
xlim([-1 3])
xticks(-1:0.5:3)
xticklabels({'-1.0','','0.0','','1.0','','2.0','','3.0'})
ylabel('$x$ [m]','interpreter','latex','FontSize',40)
ylim([-0.6 2.1])
yticks(-0.6:0.6:1.8)
yticklabels({'-0.6','0.0','0.6','1.2','1.8'})
set(gca,'TickDir','both')
xtickangle(45)



%colorbar
c = colorbar;
clim([-0.75 0.75])
c.Ticks = -0.5:0.25:0.5;
c.TickLabels = {'-0.50','-0.25','0.00','0.25','0.50'};
%c.Label.String = '$u$ [m/s]';
%c.Label.FontSize = 40;
%c.Label.Interpreter = 'latex';
%set(c.Label,'position',[4.1938-2.5 0+0.95 0],'rotation',0);
set(gca,'DataAspectRatio',[2 1 1])
crameri -roma


%print(strcat('./save_plots/velocity_',Case,'.eps'), '-depsc', '-painters');
exportgraphics(gcf, strcat('./save_plots_data/velocity_',Case,'.pdf'), 'ContentType', 'vector');
%savefig(strcat('./save_plots/velocity_',Case,'.fig'))

%set(findall(gcf,'Type','patch'),'FaceAlpha',1)
%set(gcf,'Renderer','painters')
%exportgraphics(gcf,strcat('./save_plots/velocity_',Case,'.svg'),'ContentType','vector')



% %
% figure(2)
% c = colorbar;
% clim([-0.75 0.75])
% c.Ticks = -0.5:0.25:0.5;
% c.TickLabels = {};
% %c.Label.String = '$u$ [m/s]';
% c.Label.FontSize = 40;
% c.Label.Interpreter = 'latex';
% set(c.Label,'position',[4.1938-2.5 0+0.95 0],'rotation',0);
% set(gca,'DataAspectRatio',[2 1 1])
% crameri -roma
% exportgraphics(gcf,strcat('./save_plots/velocity_colorbar.png'),'Resolution',300)



%% Plot Depth

%Go through depth and if above the shoreline make it nan

%p = polyfit(t_shore(~isnan(x_shore)),x_shore(~isnan(x_shore)),2);

for i = 1:size(Td,2)
    tm = Td(1,i);
    tav = t_shore(t_shore>tm-0.0333/2 & t_shore<tm+0.0333/2);
    xav = x_shore(t_shore>tm-0.0333/2 & t_shore<tm+0.0333/2);
    xst = mean(xav); %average shoreline over bin
    
    %xst = polyval(p,tm);


    pos_rm = Xd(:,1)>xst;
    Dd(pos_rm,i)=nan;

end

figure(1)
clf
hold on
contourf(Td-2,Xd,Dd,50,'linecolor','none')
contourf(Td,Xd,Dd,50,'linecolor','none')
contourf(Td+2,Xd,Dd,50,'linecolor','none')
plot(t_shore-2,x_shore,'k.-','linewidth',2,'markersize',10)
plot(t_shore,x_shore,'k.-','linewidth',2,'markersize',10)
plot(t_shore+2,x_shore,'k.-','linewidth',2,'markersize',10)
set(gca,'FontSize',25)

%labels and ticks
xlabel('$t$ [s]','interpreter','latex','FontSize',40)
xlim([-1 3])
xticks(-1:0.5:3)
xticklabels({'-1.0','','0.0','','1.0','','2.0','','3.0'})
ylabel('$x$ [m]','interpreter','latex','FontSize',40)
ylim([-0.6 2.1])
yticks(-0.6:0.6:1.8)
yticklabels({'-0.6','0.0','0.6','1.2','1.8'})
set(gca,'TickDir','both')
xtickangle(45)



%colorbar
c = colorbar;
clim([0 0.1])
c.Ticks = 0:0.02:0.1;
c.TickLabels = {'0.00','0.02','0.04','0.06','0.08',''};
%c.Label.String = '$h$ [m]';
%c.Label.FontSize = 40;
%c.Label.Interpreter = 'latex';
%set(c.Label,'position',[4.3813-3 0.05+0.058 0],'rotation',0);
set(gca,'DataAspectRatio',[2 1 1])
colormap(flipud(slanCM('ice')))

%set(findall(gcf,'Type','patch'),'FaceAlpha',1)
%set(gcf,'Renderer','painters')
exportgraphics(gcf, strcat('./save_plots_data/depth_',Case,'.pdf'), 'ContentType', 'vector');
%print(strcat('./save_plots/depth_',Case,'.eps'), '-depsc', '-painters');

%exportgraphics(gcf,strcat('./save_plots/depth_',Case,'.svg'),'ContentType','vector')
% 
% figure(2)
% c = colorbar;
% clim([-0.75 0.75])
% c.Ticks = -0.5:0.25:0.5;
% c.TickLabels = {};
% %c.Label.String = '$u$ [m/s]';
% c.Label.FontSize = 40;
% c.Label.Interpreter = 'latex';
% set(c.Label,'position',[4.1938-2.5 0+0.95 0],'rotation',0);
% set(gca,'DataAspectRatio',[2 1 1])
% colormap(flipud(slanCM('ice')))
% exportgraphics(gcf,strcat('./save_plots/depth_colorbar.png'),'Resolution',300)
% 
% 

%% Plot Alpha
gsT = 9.81*0.1*2;

figure(1)
clf
hold on
contourf(T-2,X,A-gsT,100,'linecolor','none')
contourf(T,X,A,100,'linecolor','none')
contourf(T+2,X,A+gsT,100,'linecolor','none')
plot(t_shore-2,x_shore,'k.-','linewidth',2,'markersize',10)
plot(t_shore,x_shore,'k.-','linewidth',2,'markersize',10)
plot(t_shore+2,x_shore,'k.-','linewidth',2,'markersize',10)
set(gca,'FontSize',25)

%labels and ticks
xlabel('$t$ [s]','interpreter','latex','FontSize',40)
xlim([-1 3])
xticks(-1:0.5:3)
xticklabels({'-1.0','','0.0','','1.0','','2.0','','3.0'})
ylabel('$x$ [m]','interpreter','latex','FontSize',40)
ylim([-0.6 2.1])
yticks(-0.6:0.6:1.8)
yticklabels({'-0.6','0.0','0.6','1.2','1.8'})
set(gca,'TickDir','both')
xtickangle(45)


%colorbar
c = colorbar;
clim([1 3])
c.Ticks = 1:0.5:3;
c.TickLabels = {'','1.5','2.0','2.5',''};
%c.Label.String = '$\alpha$ [m/s]';
%c.Label.FontSize = 40;
%c.Label.Interpreter = 'latex';
%set(c.Label,'position',[4.3813-3 0.05+3.205 0],'rotation',0);
set(gca,'DataAspectRatio',[2 1 1])
crameri -lajolla


set(findall(gcf,'Type','patch'),'FaceAlpha',1)
set(gcf,'Renderer','painters')
%print(strcat('./save_plots/alpha_',Case,'.eps'), '-depsc', '-painters');
exportgraphics(gcf, strcat('./save_plots_data/alpha_',Case,'.pdf'), 'ContentType', 'vector');

% exportgraphics(gcf,strcat('./save_plots/alpha_',Case,'.svg'),'ContentType','vector')
% 
% 
% figure(2)
% c = colorbar;
% clim([-0.75 0.75])
% c.Ticks = -0.5:0.25:0.5;
% c.TickLabels = {};
% %c.Label.String = '$u$ [m/s]';
% c.Label.FontSize = 40;
% c.Label.Interpreter = 'latex';
% set(c.Label,'position',[4.1938-2.5 0+0.95 0],'rotation',0);
% set(gca,'DataAspectRatio',[2 1 1])
% crameri -lajolla
% exportgraphics(gcf,strcat('./save_plots/alpha_colorbar.png'),'Resolution',300)



% 
% %labels and ticks
% xlabel('$t$ [s]','interpreter','latex')
% xticks(0:0.5:4)
% xticklabels({'0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0'})
% yticks(-0.4:0.4:1.6)
% yticklabels({})
% set(gca,'TickDir','both')
% set(gca,'FontSize',20)
% ylim([-0.4 1.6])
% ylim([0 0.5])
% xlim([0 2.75])
% 
% 
% %colorbar
% c = colorbar;
% clim([1.5 2.5])
% %c.Ticks = [1.25:0.25:2.25];
% %c.TickLabels = {'1.25','1.50','1.75','2.00','2.25'};
% c.Label.String = '$\alpha$ [m/s]';
% c.Label.Interpreter = 'latex';
% set(c.Label,'position',[3.1313-1.4 1.75+0.85 0],'rotation',0);
% 
% set(gca,'DataAspectRatio',[2 1 1])
% 
% return

%% Compute alpha in constant region
val = ones(size(A)); %all points start valid

p = polyfit(t_shore(~isnan(x_shore)),x_shore(~isnan(x_shore)),2);

for i = 1:size(T,2)
    tm = T(1,i);
    xs_step = polyval(p,tm);
    val(X(:,i)>xs_step,i) = nan;
end


val(T<0.75) = nan;
val(T>2) = nan;

Aval = A.*val;
val_med_A = round(median(Aval,'all','omitnan'),3)
gsT = round(9.81*0.1*2,3)

pe = abs((val_med_A-gsT)/gsT)*100;
disp(strcat("Percent Error: ",string(round(pe,2)),"%"))




if strcmp(Case,"01")
    c = 1;
elseif strcmp(Case,"02")
    c = 2;
elseif strcmp(Case,"03")
    c = 3;
end
load alpha_measure.mat
alpha_measure(c) = val_med_A
save alpha_measure.mat alpha_measure



%% Plot Dispersion

figure(1)
clf
hold on
contourf(T-2,X,D,-2.5:0.1:2.5,'LineStyle','none')
contourf(T,X,D,-2.5:0.1:2.5,'LineStyle','none')
contourf(T+2,X,D,-2.5:0.1:2.5,'LineStyle','none')
plot(t_shore-2,x_shore,'k.-','linewidth',2,'markersize',10)
plot(t_shore,x_shore,'k.-','linewidth',2,'markersize',10)
plot(t_shore+2,x_shore,'k.-','linewidth',2,'markersize',10)
set(gca,'FontSize',25)

%labels and ticks
xlabel('$t$ [s]','interpreter','latex','FontSize',40)
xlim([-1 3])
xticks(-1:0.5:3)
xticklabels({'-1.0','','0.0','','1.0','','2.0','','3.0'})
ylabel('$x$ [m]','interpreter','latex','FontSize',40)
ylim([-0.6 2.1])
yticks(-0.6:0.6:1.8)
yticklabels({'-0.6','0.0','0.6','1.2','1.8'})
set(gca,'TickDir','both')
xtickangle(45)

%colorbar
c = colorbar;
clim([-2.5 2.5])
c.Ticks = -2:1:2;
c.TickLabels = {'-2','-1',' 0',' 1',' 2'};
%c.Label.String = '$\partial u / \partial x$ [s$^{-1}$]';
%c.Label.Interpreter = 'latex';
%c.Label.FontSize = 40;
%set(c.Label,'position',[4.3813-4.5 0.05+3 0],'rotation',0);
set(gca,'DataAspectRatio',[2 1 1])
crameri bam

set(findall(gcf,'Type','patch'),'FaceAlpha',1)
set(gcf,'Renderer','painters')
%print(strcat('./save_plots/dudx_',Case,'.eps'), '-depsc', '-painters');
exportgraphics(gcf, strcat('./save_plots_data/dudx_',Case,'.pdf'), 'ContentType', 'vector');

% exportgraphics(gcf,strcat('./save_plots/dudx_',Case,'.svg'),'ContentType','vector')
% 
% figure(2)
% c = colorbar;
% clim([-0.75 0.75])
% c.Ticks = -0.5:0.25:0.5;
% c.TickLabels = {};
% %c.Label.String = '$u$ [m/s]';
% c.Label.FontSize = 40;
% c.Label.Interpreter = 'latex';
% set(c.Label,'position',[4.1938-2.5 0+0.95 0],'rotation',0);
% set(gca,'DataAspectRatio',[2 1 1])
% crameri -tokyo
% exportgraphics(gcf,strcat('./save_plots/dudx_colorbar.png'),'Resolution',300)


%% Get swawsh tip length
clear x_tip t_tip

x = X(:,1);
t = [T(1,:) T(1,:)+2];

step = 1;

if strcmp(Case,"01")
    t_start = 0.7;
    t_end = 1.4;
elseif strcmp(Case,"02")
    t_start = 0.9;
    t_end = 1.8;
elseif strcmp(Case,"03")
    t_start = 1;
    t_end = 2.25;
end


for j = find(t>t_start & t<=t_end) %C3: (1 2.3), %C2: (0.9 1.8]
    if j > size(T(1,:),2)
        i = j - size(T(1,:),2);
        down = 1;
    else
        i = j;
        down = 0;
    end

    %D(:,i);
    figure(2)
    clf
    d = medfilt1(D(:,i));
    plot(x,d)
    hold on
    yline(0)

    ind = find(diff(sign(medfilt1(D(:,i))))==-2); %index where sign changes
    if strcmp(Case,"03")
        if i == 17
            ind = ind(2);
        elseif i == 18
             ind = find(diff(sign(medfilt1(D(:,i)-0.1)))==-2);
        elseif i == 19
            ind = ind(1);
        end
    elseif strcmp(Case,"02")
        if i == 18
            ind = find(diff(sign(medfilt1(D(:,i)-0.1)))==-2);
        elseif i == 19
            ind = find(diff(sign(medfilt1(D(:,i)-0.15)))==-2);
        end
    elseif strcmp(Case,"01")
        if i == 10
            ind = find(diff(sign(medfilt1(D(:,i)-0.15)))==-2);
        elseif i == 11
            ind = find(diff(sign(medfilt1(D(:,i)-0.1)))==-2);
        elseif i == 12
            ind = ind(2);
        end
    end
    
    if ~isempty(ind)
        if strcmp(Case,"03")
            if i == 18
                 x_tip(step) = interp1([d(ind)-0.1 d(ind + 1)-0.1],[x(ind) x(ind + 1)],0);
            else
                x_tip(step) = interp1([d(ind) d(ind + 1)],[x(ind) x(ind + 1)],0);
            end
        elseif strcmp(Case,"02")
            if i == 18
                 x_tip(step) = interp1([d(ind)-0.1 d(ind + 1)-0.1],[x(ind) x(ind + 1)],0);
            elseif i == 19
                 x_tip(step) = interp1([d(ind)-0.15 d(ind + 1)-0.15],[x(ind) x(ind + 1)],0);
            else
                x_tip(step) = interp1([d(ind) d(ind + 1)],[x(ind) x(ind + 1)],0);
            end
        elseif strcmp(Case,"01")
            if i == 10
                x_tip(step) = interp1([d(ind)-0.15 d(ind + 1)-0.15],[x(ind) x(ind + 1)],0);
            elseif i == 11
                x_tip(step) = interp1([d(ind)-0.1 d(ind + 1)-0.1],[x(ind) x(ind + 1)],0);
            else
                x_tip(step) = interp1([d(ind) d(ind + 1)],[x(ind) x(ind + 1)],0);
            end
        end
        plot(x_tip(step),0,'ro','linewidth',2)
    else
        x_tip(step) = nan;
    end

    %title(strcat(string(j),"/",string(max(find(t>1 & t<2.2)))))
    drawnow()



    if down == 0
        t_tip(step) = t(i);
    elseif down == 1
        t_tip(step) = t(i) + 2;
    end

    step = step + 1;

end

%
figure(1)
hold on
plot(t_tip,x_tip,'k:','linewidth',2)

exportgraphics(gcf, strcat('./save_plots_data/dudx_',Case,'_tip.pdf'), 'ContentType', 'vector');

%% tip width
x_tip;
t_tip;

x_shore;
t_shore;


for i = 1:length(t_tip)
    ind = find(abs(t_shore-t_tip(i)) == min(abs(t_shore-t_tip(i))),1); %closest shore index
    w(i) = x_shore(ind) - x_tip(i);
end

max(w); %max tip width
max(x_shore); %max shoreline location

if strcmp(Case,"01")
    c = 1;
elseif strcmp(Case,"02")
    c = 2;
elseif strcmp(Case,"03")
    c = 3;
end
load tip_width.mat
tip_width_dimless(c) = max(w)/max(x_shore)
save tip_width.mat tip_width_dimless

load max_shoreline.mat
xs_max(c) = max(x_shore)
ts_max(c) = t_shore(find(x_shore == max(x_shore),1))
save max_shoreline.mat xs_max ts_max