%==========================================================================
%==========================================================================

clear;
clc;

%% Setup
% All for wave angle = 0 degrees
Case = '03'; %Define wave case
path_comb = strcat("R:\Ben Davidson\Queens 2024 Data\S_A00_C123\00",Case,"_comb_S");
path_cc = "R:\Ben Davidson\Queens 2024 Code\Image Preprocessing";
path_swl = "R:\Ben Davidson\Queens 2024 Data\background";


%% Maybe add file check?

%% LOAD
load WG_locs.mat %WG LOCATIONS
%load files combined
load(strcat(path_comb,'\',"final_shoreline_ea_comb.mat"),'x_shore','t_shore');
%Load files both
load(strcat(path_cc,'\cc_smooth.mat'))
load(strcat(path_swl,'\A00_SWL.mat'))
x_ref = X_swl;

load(strcat(path_comb,"\combined_wave_gauges.mat"))
%% Processes
% Need 4 x locations of sensors
x = [A5_px(1) B5_px(1) B6_px(1) A6_px(1)];


max_xs = max(x_shore);
% Add x location of max run-up as last point
x = [x max_xs];

%convert x locations to distance in meters relative to SWL
x = (x_ref-x)*mm_px/1000;

%save depth over space
xd = -0.6:0.05:1.8;
td = linspace(0,2,61);

[T,X] = meshgrid(td,xd);
D = nan(size(T));


%figure(1)
%clf

figure(2)
clf
set(gca, 'Units','normalized', 'Position',[0.15 0.2 0.75 0.7]);
tiledlayout(2,3,"TileSpacing","compact")
index = 1;

for i = 1:61
    %figure(1)
    %subplot(2,3,index)
    %clf
    %hold on
    %plot([-0.6 1.8],[-0.6 1.8]*(1/10),'k-','linewidth',2)
    %plot([-0.6 0],[0 0],'b--','linewidth',1)
    %ylim([-0.06 0.6])
    %xlim([-0.6 1.8])
    %xlabel('x [m]')
    %ylabel('z [m]')
    %set(gca,'FontSize',20)
    
    % %sensors
    % for se = 1:length(x)
    %     xline(x(se),'r--')
    % end

    d = [deA{2}(i,2) deB{5}(i,2) deB{6}(i,2) deA{3}(i,2) 0]; %[A5 B5 B6 A6]
    y = d + x/10;
    
    %plot(x(1),y(1),'ro','markerfacecolor','r')

    %plot(x(2),y(2),'ro','markerfacecolor','r')
    %plot(x(3),y(3),'ro','markerfacecolor','r')

    %plot(x(4),y(4),'ro','markerfacecolor','r')
    %plot(x(5),y(5),'ro','markerfacecolor','r')



    % Model function: y = a*x.^2 + b*x + c
    modelFun = @(b,x) b(1)*x.^2 + b(2)*x + b(3);
    
    % Initial guess for parameters [a, b, c]
    beta0 = [1, 1, 1];
    
    % Fit
    [beta, ~, ~, ~, MSE, ~] = nlinfit(x, d, modelFun, beta0);

    err = sqrt(MSE);

    %plot fit
    xfit = -0.6:0.1:1.8;
    dfit = modelFun(beta, xfit); %local depth

    %plot_err(xfit,dfit,err,max(x))
    % % % %dfit(dfit<0) = 0; %set negative depth as zero
    % % % yfit = dfit + xfit/10; %fit in y coords.
    % % % plot(xfit, yfit, 'r-')
    % % % 
    % % % ME = sqrt(MSE); %mean error -> av. dist. between fit and actual points
    % % % 
    % % % %plot errorbars
    % % % xer = -0.5:0.1:1.5;
    % % % yer = modelFun(beta,xer)+xer/10;
    % % % errorbar(xer,yer,errors,'o')

    %drawnow()

    %Save Depth Fit
    dd = modelFun(beta,xd);
    dd(dd<0) = nan;
    dd(xd>max(x)) = nan;
    D(:,i) = dd;

    if any(i == [1 13 25 37 49 61])
        nexttile
        hold on
        plot([-0.6 1.8],[-0.6 1.8]*(1/10),'k-','linewidth',2)
        plot([-0.6 0],[0 0],'b--','linewidth',1)
        ylim([-0.06 0.5])
        xlim([-0.6 1.8])
        xticks([])
        yticks([])

        %plot sensors
        plot(x(1),y(1),'ro','markerfacecolor','r')
    
        plot(x(2),y(2),'ro','markerfacecolor','r')
        plot(x(3),y(3),'ro','markerfacecolor','r')
    
        plot(x(4),y(4),'ro','markerfacecolor','r')
        plot(x(5),y(5),'ro','markerfacecolor','r')

        %plot depth
        plot_err(xfit,dfit,err,max(x))

        if any(index == [1 4]) %left side
            ylabel('$z$ [m]', 'interpreter','latex')
            yticks([0:0.1:0.5])
            
        end
        if any(index == [4 5 6]) %bottom
            xlabel('$x$ [m]', 'interpreter','latex')
            xticks([-0.5:0.5:1.5])
        end
        
       
        set(gca,'FontSize',25)
        set(gca,'TickLabelInterpreter','Latex')
        title(strcat("$t$ = ",string(td(i)), " s"),'interpreter','latex')
    
        index = index + 1;
    end



end

return

%% fit depth over space
figure(2)
clf
hold on
contourf(T,X,D,50,'linecolor','none')
plot(t_shore,x_shore,'k-','linewidth',2)

colorbar
caxis([0 0.15])
colormap sky

xlabel('t [s]')
ylabel('x [m]')
set(gca,'FontSize',20)

%add shoreline
plot(t_shore,x_shore,'k-.','linewidth',2)


Td = T;
Xd = X;
Dd = D;


for i = 1:size(D,2)
    ts = T(1,i);
    xs = mean(x_shore(t_shore>ts-0.0333/2 & t_shore<ts+0.0333/2));
    D(X(:,i)>xs,i) = nan;
end

figure(3)
clf
contourf(T,X,D,50,'linecolor','none')
hold on
plot(t_shore,x_shore,'k-','linewidth',2)

colorbar
caxis([0 0.15])
colormap sky

xlabel('t [s]')
ylabel('x [m]')
set(gca,'FontSize',20)


Td = T;
Xd = X;
Dd = D;

%%

save(strcat(path_comb,'\ea_depth_field_comb.mat'),'Td','Xd','Dd')




function plot_err(xfit,dfit,err,xmax)
    dfit(dfit<0) = 0;


    d_high = dfit + err;
    d_low = dfit - err;

    d_high(d_high<0) = 0;
    d_low(d_low<0) = 0;

    dfit(xfit>xmax) = 0;
    d_high(xfit>xmax) = 0;
    d_low(xfit>xmax) = 0;

    y_high = d_high + xfit/10;
    y_low = d_low + xfit/10;
    yfit = dfit + xfit/10;

    yp = [y_high fliplr(y_low)];
    xp = [xfit fliplr(xfit)];

    %plot main
    plot(xfit,yfit,'b-','linewidth',1)

    %plot error
    fill(xp,yp,'r','facecolor','cyan','facealpha',0.25,'edgecolor','none')

end
