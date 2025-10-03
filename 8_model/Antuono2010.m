%==========================================================================
% Antuono2010.m
%
% Solve Antuono 2010 (non-dimensional).
% 
% Author: B. Davidson
% Last Updated: 3 October 2025
%==========================================================================

clear;
close all;
clc;

%% Setup
alpha2 = 2.4607; %alpha sea-ward of the shock (non-dimensional)
%alpha1 = 2; %alpha shore-ward of the shock


%Computation Grid
t_end = 5; %simulation end time
% 
% dt = 0.005;
% dx = 0.001;
% 
% tm = 0:dt:t_end;
% xm = -1:dx:2;

%Runge-Kutta timestep details
tau0 = 10^(-4);
n = 12;

%% Shock Solution - Bore approaching SWL

%HY Code
% [t_s_mat, x_s_mat, beta2_mat, u2_mat, s_mat, d2_mat, t_collapse, index_collapse, alpha2_mat] = function(alpha2,alpha1,tau0,n,t_end);
% Us = u2_mat(index_collapse);

% INITIALIZE
ti = 1; %time index starts at one
t_mat(ti) = 0; %initial time is zero
x_s_mat(ti) = -1; %initialize the bore at x = -1, the toe of the beach

% water is still before the shock enters sloping beach: u1 = 0 and d1 = sqrt(-x)
%       at x = -1; d1 = 1; c1 = 1; since d1 = sqrt(c1)
d1 = -x_s_mat(ti); %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check with HS
d1_mat(ti) = d1;
c1 = sqrt(d1);          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check with HS

u1 = 0; %cross-shore velocity - shoreward of the shock wave -- must be zero since water is still before the shock (velocity infront of shock)
u1_mat(ti) = u1;

t1 = t_mat(ti);

alpha1 = u1 + 2*c1 + t1;

z = eq_5p3_solve(c1,alpha1,alpha2); %solve Antuono Equation(5.3) to get d2
d2_mat(ti) = z^2*c1;


d2 = d2_mat(ti);
s(ti) = u1 + sqrt(0.5*(d2^2/d1 + d2)); %shock wave velocity (Eq. 3.9)

u2_mat(ti) = alpha2 - 2*sqrt(d2_mat(ti)) - t_mat(ti); %from def'n of alpha (between 5.2 and 5.3)
beta2_mat(ti) = u2_mat(ti) -2*sqrt(d2_mat(ti)) + t_mat(ti); %(Eq. 2.2)

while d1>10^(-12)
    %numerical integration of Antuono 2010 Equation (5.2) with RK-4
    dt_rk4 = tau0*(1 - (1 - d1)^n); %RK-4 step size

    %% RK-4 Integration

    % initial x, t...
    x_n0 = x_s_mat(ti);
    t_n0 = t_mat(ti);
    s_n0 = s(ti);

    %First Step
    k1 = s_n0; %s evaluated at x_n0 and t_n0

    %Second Step
    t_step2 = t_n0 + dt_rk4/2;
    x_step2 = x_n0 + dt_rk4*(k1/2);

    d1 = -x_step2;      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check with HS
    c1 = sqrt(d1);
    alpha1 = u1 + 2*c1 + t_step2;
    z = eq_5p3_solve(c1,alpha1,alpha2);
    d2 = z^2*c1;
    k2 = u1 + sqrt(0.5*(d2^2/d1 + d2));

    %Third Step
    t_step3 = t_step2;
    x_step3 = x_n0 + dt_rk4*(k2/2);

    d1 = -x_step3;       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check with HS
    c1 = sqrt(d1);
    alpha1 = u1 + 2*c1 + t_step3;
    z = eq_5p3_solve(c1,alpha1,alpha2);
    d2 = z^2*c1;
    k3 = u1 + sqrt(0.5*(d2^2/d1 + d2));

    %Fourth Step
    t_step4 = t_n0 + dt_rk4;
    x_step4 = x_n0 + dt_rk4*k3;
    
    d1 = -x_step4;      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check with HS
    c1 = sqrt(d1);
    alpha1 = u1 + 2*c1 + t_step4;
    z = eq_5p3_solve(c1,alpha1,alpha2);
    d2 = z^2*c1;
    k4 = u1 + sqrt(0.5*(d2^2/d1 + d2));
    

    x_n1 = x_n0 + (dt_rk4/6)*(k1 + 2*k2 + 2*k3 + k4);
    t_n1 = t_n0 + dt_rk4;

    %% Update timestep from RK-4 Integration
    ti = ti + 1;
    t_mat(ti) = t_n1;
    x_s_mat(ti) = x_n1;
    d1 = -x_n1;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check with HS
    d1_mat(ti) = d1;
    c1 = sqrt(d1);
    alpha1 = u1 + 2*c1 + t_n1;
    z = eq_5p3_solve(c1,alpha1,alpha2);
    d2 = z^2*c1;
    d2_mat(ti) = d2;
    s(ti) = u1 + sqrt(0.5*(d2^2/d1 + d2));

    u2_mat(ti) = s(ti) - sqrt(0.5*((d1^2/d2)+d1));
    u1_mat(ti) = u1;
    beta2_mat(ti) = u2_mat(ti) -2*sqrt(d2_mat(ti)) + t_mat(ti);
    
end

alpha2_mat = t_mat + 2*sqrt(d2_mat) + u2_mat;

xs0 = x_s_mat;%still shoreline
xs0(x_s_mat<0) = 0;
xs0(x_s_mat>0) = nan;

%% Swashzone

dt_swash = tau0; %timestep to solve the swash solution
t0 = x_s_mat(end)*((t_mat(end-1)-t_mat(end))/(x_s_mat(end)-x_s_mat(end-1)))+t_mat(end); %time shock reaches shoreline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Why is mine diff. from HYS?
x0 = 0; %position when shock reaches shoreline - should maybe be zero?
% collapse_ind = length(t_mat)+1;
% t_mat = [t_mat t0];
% x_s_mat = [x_s_mat x0];
% s(end+1) = alpha2-t0;

t_swash = [t0 t0+dt_swash:dt_swash:t_end];
xs_swash = x0 + alpha2*(t_swash-t0)-(t_swash.^2)/2+(t0^2)/2; %Eq. 6.2
s_swash = alpha2-t0;
u2_swash = alpha2-t_swash;
alpha2_swash = u2_swash + t_swash; % At the shoreline, C = 0;

xb_a = [x_s_mat xs_swash];
xs_a = [xs0 xs_swash];
ts_a = [t_mat t_swash];
us_a = [u2_mat u2_swash];
alpha2_a = [alpha2_mat alpha2_swash];
beta2_swash = u2_swash + t_swash; % At the shoreline, C = 0;
beta2_a = [beta2_mat beta2_swash];

zero_index = find(ts_a == t0);
t0;
Us = us_a(zero_index);

%% Solving behind shock


t_vect = 0:0.01:4;
ind = 1;

i_plotsB = [1000:1000:25000 30000];

for i = 1:length(i_plotsB)
    ind = i_plotsB(i);
    %(ts,xs) ordered pair of beta curve intersecting bore
    ts = ts_a(ind); %bore time inquery
    xs = xb_a(ind); %bore space inquery

    Beta_pairs = beta2_a(ind);

    %at each bore time/space coordinate (ts,xs) at any time, we can solve
    %for the position of the point along the constant Beta curve
    x_vect_B(:,i) = xs + ((alpha2 + 3*Beta_pairs)/4).*(t_vect-ts) - (t_vect.^2)/2 + (ts^2)/2;

    %x_vect_A(:,i) = xs + ((3*alpha2 + Beta_pairs)/4).*(t_vect-ts) - (t_vect.^2)/2 + (ts^2)/2;

    %valid = x_vect(:,i)<xs;
    validB = t_vect>ts;
    %validA = t_vect<ts;
    %valid(t_vect<ts) = 0;

    valB(:,i) = validB;
    %valA(:,i) = validA;

end

% alpha does not start along the bore, rather starts at the offshore
% boundary

t_vect; %same
x_start = -1;

t_start = 0.5; %start time inquery
t_vectA = 0.5:0.1:4;
for i = 1:length(t_start)
    x_vect_A(:,i) = x_start + ((3*alpha2 + Beta_pairs)/4).*(t_vectA-t_start(i)) - (t_vectA.^2)/2 + (t_start(i)^2)/2;
end



%% Once we have Beta at any point, can integrate to get alpha curves...


dt = 0.005;
dx = 0.001;

tm = 0:dt:t_end;
xm = -1:dx:2;

U = nan(length(tm),length(xm));
D = nan(length(tm),length(xm));
B = nan(length(tm),length(xm));

[X,T] = meshgrid(xm,tm);

%loop through time on my grid
for ti = 1:length(tm)
    tt = tm(ti); %time at step
    x_shock = interp1(ts_a,xb_a,tt); %shock position at this time step

    %calculate x positions - and list with corresponding beta values for
    %points on characteristic line
    x_char = xb_a + ((alpha2 + 3*beta2_a)/4).*(tt-ts_a) - (tt^2)/2 + (ts_a.^2)/2;
    b_char = beta2_a;

    %trim x_char and b_char to only be the points behind x_shock
    i_trim = find(x_char>x_shock,1,'first'); %find first point where x is larger than shock location
    i_trim = i_trim - 1; %go one step back to be only behind shock

    x_char = x_char(1:i_trim);
    b_char = b_char(1:i_trim);


    valid = xm<x_shock;

    if sum(valid)==0 & length(x_char)<1 %if no valid points and no points in x_char - SKIP
        continue
    end



    b_interp = interp1(x_char,b_char,xm(valid));

    B(ti,valid) = b_interp;


end

%% Compute U and D

U = (alpha2 + B)/2 - T;
D = ((alpha2-B)/4).^2;
Eta = D + X;

%if x is behind 0, then if U and D are nan, then the values are zero
%%%%%%%%%%% NEED TO EDIT: this is only true at the start not the end...
% for i = 1:length(tm)
%     Ustp = U(i,:);
%     Ustp(isnan(Ustp) & xm<=0) = 0;
%     U(i,:) = Ustp;
% 
% 
%     Dstp = D(i,:);
%     Dstp(isnan(Dstp) & xm<=0) = 0;
%     D(i,:) = Dstp;
% end
%% Package up Antuono model to run over experiment data
    % Want x, t, u, d, maybe eta, alpha, Us, t0, x0 (always zero)

    
save(strcat('alpha2_',string(alpha2),'.mat'))
disp('DONE')
return


%%
l0 = 1; %[m]
g = 9.81;
s = 1/10;
t0d = sqrt(l0/(g*s));
u0 = sqrt(g*s*l0);


figure(1)
clf
subplot(2,2,4)
ind = find(round(xm,4)==-1); %0.01
plot((tm-t0)*t0d,U(:,ind)*u0,'k-')

xlabel('t')
ylabel('u')
set(gca,'FontSize',15)
% xlim([0 2])
% xticks(0:0.5:2)
% ylim([-1 1])
% yticks(-1:0.5:1)


subplot(2,2,2)
plot((tm-t0)*t0d,D(:,ind)*(s)*(l0),'k-')
title('B')
xlabel('t')
ylabel('d')
set(gca,'FontSize',15)
%xlim([0 2])
%xticks(0:0.5:2)
% ylim([0 0.4])
% yticks(0:0.1:0.4)


subplot(2,2,3)
ind = find(round(xm,4)==0.5);
plot((tm-t0)*t0,U(:,ind)*1.5,'k-')

xlabel('t')
ylabel('u')
set(gca,'FontSize',15)
xlim([0 2])
xticks(0:0.5:2)
ylim([-1 1])
yticks(-1:0.5:1)


subplot(2,2,1)
plot((tm-t0)*t0,D(:,ind)*0.65,'k-')
title('A')
xlabel('t')
ylabel('d')
set(gca,'FontSize',15)
xlim([0 2])
xticks(0:0.5:2)
ylim([0 0.4])
yticks(0:0.1:0.4)


%%
l0 = 1; %[m]
g = 9.81;
s = 1/10;
t0d = sqrt(l0/(g*s));
u0 = sqrt(g*s*l0);


figure(1)
clf
subplot(2,1,2)
ind = find(round(xm,4)==-0.223); %0.01
plot((tm),U(:,ind),'k-')

xlabel('t')
ylabel('u')
set(gca,'FontSize',15)
% xlim([0 2])
% xticks(0:0.5:2)
% ylim([-1 1])
% yticks(-1:0.5:1)


subplot(2,1,1)
plot((tm),D(:,ind),'k-')
title('x = -1')
xlabel('t')
ylabel('d')
set(gca,'FontSize',15)
%xlim([0 2])
%xticks(0:0.5:2)
% ylim([0 0.4])
% yticks(0:0.1:0.4)

%%
figure(3)
clf

for i = 1:length(tm)
    plot([-1 2],[-1 2],'linewidth',2)
    hold on
    plot([-1 0],[0 0],'k:','linewidth',1)
    %plot(xm,D(i,:),'k-')
    plot(xm,Eta(i,:),'k-','linewidth',2)
    xlim([min(xm) max(xm)])
    ylim([-1 2])
    xlabel('x [-]')
    ylabel('y [-]')
    set(gca,'FontSize',15)

    drawnow();
    if round(tm(i),4) == round(t0,2)
        %return
    end
    %pause(0.1)
    hold off
end


%%

%% Shoreline

g = 9.81;
s = 1/10;


figure(2)
clf
tp = (ts_a-t0)*2.57;
xp = (xs_a)*6.5;
plot(tp,xp,'-','linewidth',2)
xlabel('t [s]')
ylabel('x [m]')
set(gca,'FontSize',15)
ylim([0 2])
xlim([0 4])




%%


function z = eq_5p3_solve(c1,alpha1,alpha2)
    %set up Equation 5.3 (Antuono 2010) as polynomial
    z6 = 1;
    z5 = 0;
    z4 = -9*c1;
    z3 = 8*sqrt(c1)*(alpha2 - alpha1 + 2*c1);
    z2 = -1*(2*(alpha2 - alpha1 + 2*c1)^2+c1^2);
    z1 = 0;
    z0 = c1^3;
    p = [z6 z5 z4 z3 z2 z1 z0];

    %solve for roots of z
    r = roots(p);

    %the root must be real
    r_real = r(imag(r)==0);
    %"the correct root is that satisfying z>=sqrt(c1)
    z = r_real(r_real>=sqrt(c1));

    %error if no solutions or more than 1 solution
    if isempty(z) | length(z)>1
        error('error')
    end
end


   