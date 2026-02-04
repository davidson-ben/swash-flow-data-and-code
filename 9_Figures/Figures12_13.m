clear;
%close all;
clc;

%plot dudx at specific cross-shore locations
Case = "03";
path_comb = strcat("R:\Ben Davidson\Queens 2024 Data\S_A00_C123\00",Case,"_comb_S");

%load(strcat(path_comb,'/ea_depth_field.mat'),'Xd','Td','Dd')
load(strcat(path_comb,'/plot_uhad.mat'))

%% Plot Dispersion

p = polyfit(t_shore(~isnan(x_shore)),x_shore(~isnan(x_shore)),2);
figure(1)
clf
hold on
plot(t_shore,x_shore,'b.-')
tf = 0:0.1:2.5;
xf = polyval(p,tf);
plot(tf,xf,'k-')

%%
D0 = nan(size(D));
D1 = D;
D2 = nan(size(D));
for i = 1:size(D1,2) %step through time
    ct = T(1,i); %current time
    if ct < 1.2 %anything above the shoreline moves to the current cycle under the shoreline
        xs = polyval(p,ct); %current shoreline position
        D2(X(:,i)>xs,i) = D1(X(:,i)>xs,i);
        D1(X(:,i)>xs,i) = nan;
    end

    first_shoreline = find(~isnan(x_shore),1); %first non-nan shoreline position
    if ct > t_shore(first_shoreline)+2 %anything under the 2nd shoreline moves under the first shoreline
        xs = polyval(p,ct-2); %current shoreline position
        D0(X(:,i)<xs,i) = D1(X(:,i)<xs,i);
        D1(X(:,i)<xs,i) = nan;
    end
end

Tc = [T-2 T T+2];
Dc = [D0 D1 D2];
Xc = [X X X];

figure(1)
clf
hold on
contourf(Tc,Xc,Dc,-2.5:0.1:2.5,'LineStyle','none')
plot(t_shore,x_shore,'k.-','linewidth',2,'markersize',10)
plot(t_shore+2,x_shore,'k.-','linewidth',2,'markersize',10)
set(gca,'FontSize',25)

%labels and ticks
xlabel('$t$ [s]','interpreter','latex','FontSize',40)
xlim([-1 3])
xticks(-1:1:3)
xticklabels({'-1.0','0.0','1.0','2.0','3.0'})
ylabel('$x$ [m]','interpreter','latex','FontSize',40)
ylim([-0.6 1.8])
yticks(-0.6:0.6:1.8)
yticklabels({'-0.6','0.0','0.6','1.2','1.8'})
set(gca,'TickDir','both')

%colorbar
c = colorbar;
clim([-2.5 2.5])
c.Ticks = -2:1:2;
c.TickLabels = {'-2','-1',' 0',' 1',' 2'};
c.Label.String = '$\partial u / \partial x$ [s$^{-1}$]';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 40;
set(c.Label,'position',[4.3813-4.5 0.05+3 0],'rotation',0);
set(gca,'DataAspectRatio',[2 1 1])
%crameri -tokyo

%set(findall(gcf,'Type','patch'),'FaceAlpha',1)
%set(gcf,'Renderer','painters')
%exportgraphics(gcf,strcat('./save_plots/dudx_',Case,'.svg'),'ContentType','vector')


%%
figure(2)
clf
hold on

%xcs = Xc(:,1);
%x = 0.6; %[-0.75:0.05:1.50]
%ind = find(abs(xcs-x)==0);

mt = 0:0.01:3;
md = 2./(3*mt);

xcs = Xc(:,1);
col = parula(length(xcs));
fil = 1; %filter
n = 3;
%w = nan(size(xcs));
%w(xcs<0) = 1;
%w(xcs>=0) = 2;
%w(xcs>=0.6) = 1;

ll = -0.0; %lower limit
ul = 0.6; %upper limit

for ind = 1:46
    if fil == 1
        if xcs(ind)<ll
            %plot(Tc(ind,:),medfilt1(Dc(ind,:),n),'-','color',col(ind,:),'LineWidth',1)
        elseif xcs(ind)>=ll && xcs(ind)<ul
            plot(Tc(ind,:),medfilt1(Dc(ind,:),n),'-','color',col(ind,:))
        else
            %plot(Tc(ind,:),medfilt1(Dc(ind,:),n),'-','color',col(ind,:),'LineWidth',1)
        end
    elseif fil == 0
        if xcs(ind)<ll
            plot(Tc(ind,:),Dc(ind,:),'-','color',col(ind,:),'LineWidth',1)
        elseif xcs(ind)>=ll && xcs(ind)<ul
            plot(Tc(ind,:),Dc(ind,:),'.-','color',col(ind,:),'LineWidth',2,'markeredgecolor','k')
        else
            plot(Tc(ind,:),Dc(ind,:),'-','color',col(ind,:),'LineWidth',1)
        end
    end
end
plot(mt,md,'k--')
set(gca,'YScale','log')
set(gca,'XScale','log')

c = colorbar;
clim([min(xcs) max(xcs)])

xlabel('time [s]')
ylabel('du/dx [s^{-1}]')
set(gca,'FontSize',20)

c.Label.String = "Cross-shore position [m]";
c.Ticks = -0.75:0.25:1.5;
c.TickLabels = {'-0.75','-0.50','-0.25','0.00','0.25','0.50','0.75','1.00','1.25','1.50'};



%exportgraphics(gcf,strcat("R:\Ben Davidson\Figures\dudx_lines_",Case,'.jpg'))


%%
figure(3)
clf
hold on

%xcs = Xc(:,1);
%x = 0.6; %[-0.75:0.05:1.50]
%ind = find(abs(xcs-x)==0);

mt = 0:0.01:3;
md = 2./(3*mt);

xcs = Xc(:,1);
col = parula(length(xcs));
fil = 1; %filter
n = 3;
%w = nan(size(xcs));
%w(xcs<0) = 1;
%w(xcs>=0) = 2;
%w(xcs>=0.6) = 1;

ll = 0.0; %lower limit
ul = 0.6; %upper limit

for ind = 1:46
    if fil == 1
        if xcs(ind)<ll
            %plot(Tc(ind,:),medfilt1(Dc(ind,:),n),'-','color',col(ind,:),'LineWidth',1)
        elseif xcs(ind)>=ll && xcs(ind)<ul
            plot(Tc(ind,:),medfilt1(Dc(ind,:),n),'-r','LineWidth',0.5)
        else
            %plot(Tc(ind,:),medfilt1(Dc(ind,:),n),'-','color',col(ind,:),'LineWidth',1)
        end
    elseif fil == 0
        if xcs(ind)<ll
            plot(Tc(ind,:),Dc(ind,:),'-','color',col(ind,:),'LineWidth',1)
        elseif xcs(ind)>=ll && xcs(ind)<ul
            plot(Tc(ind,:),Dc(ind,:),'.-','color',col(ind,:),'LineWidth',2,'markeredgecolor','k')
        else
            plot(Tc(ind,:),Dc(ind,:),'-','color',col(ind,:),'LineWidth',1)
        end
    end
end
plot(mt,md,'k--','linewidth',3)
set(gca,'YScale','log')
set(gca,'XScale','log')
xlim([0.1 3])
ylim([0.001 10])



xlabel('time [s]','interpreter','latex')
ylabel('$\partial u/ \partial x [s^{-1}]$','interpreter','latex')
set(gca,'FontSize',30)
set(gca,'TickLabelInterpreter','latex')

exportgraphics(gcf,strcat("R:\Ben Davidson\Figures\dudx_lines_",Case,'.pdf'),'ContentType','vector')
