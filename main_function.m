%main function for tracking FRF and higher dimension curves based on the
%system and the system parameters

clear;
clc;
close all;

global Fw alpha zeta
Fw=0.05;
alpha=1; %nonlinear-linear ratio ksn*F0^2/ks^3, can be studied to see the impact of the 
%nonlinear of the primary system on the NERS-TMDI performance
zeta=0.01; %damping in the primary system

ga=7.4704e-04;%3.7037e-04;%7.4704e-04

%system parameters structure sysP=[b mu lambda gamma sigma kappa rho]
sysP=[0.0423 0.04 1*0.0596 ga 1*0.379 1*1.5+0*1.0161 1*0.0618 0.000];

%FRF construct the resonant curve (main) and detect fold bifurcation
%the user can decide to track another fold bifurcation or stop and move
%to the next continuation path
x_mrc=FRF(sysP); 
x_LPC=x_mrc(:,end-1);%Limit Point Cycle
x=x_mrc;
mx1=(sqrt(x(1,:).^2+x(2,:).^2+x(3,:).^2+x(4,:).^2+x(5,:).^2)); %primary sys disp. amplitude
x_m1=islocalmax(mx1); %detect local maximas of the mrc
idx=find(x_m1); %indexes for maximas of the mrc
for hh=1:length(idx)
    x_max(hh,:)=[x_mrc(:,idx(hh))',x(end,idx(hh)),mx1(idx(hh))]; %maximum for the MRC
end

%compute bifurcation diagrams continuing maximas of MRC
prompt="compute bifurcation diag.? 'y' for yes, 'n' for no: ";
enter=input(prompt,'s');
if enter=='y'
    enter=true;
    pl0=true;
else
    enter=false;
    pl0=false;
end
while enter
F_init=Fw;


if length(idx)==1
    x1=Bif(x_max(1,:),sysP,F_init);
elseif length(idx)==2
    x1=Bif(x_max(1,:),sysP,F_init);
    x2=Bif(x_max(2,:),sysP,F_init);
elseif length(idx)==3
    x1=Bif(x_max(1,:),sysP,F_init);
    x2=Bif(x_max(2,:),sysP,F_init);
    x3=Bif(x_max(3,:),sysP,F_init);
end
enter=false;
end
%L1 represents the first level of continuation where the Limit Point Cycle
%detected at the fold bif. It is continued based on the user's input. For
%this analysis for Level 1, we track the Limit Point as a function of
%forcing Fw.

%% ----------------------Level 1------------------------------------
prompt="do you want to conduct a level 1 continuation? ";
q=input(prompt,'s');
if q=='y'
    %2 direction, user choose a positive or negative step to impose a direction
%to the continuation path
fprintf('\n This is the continuation for level 1 \n\n')
fprintf(' ')
prompt="Choose a positive or negative step value: ";
step=input(prompt);
%direction 1
x_L1_d1 = L1(x_LPC,sysP,step);
%x_dF=x_L1_d1(:,end);%Turning point
xt=x_L1_d1;
xt=[xt(:,1) xt(:,3:end)];
%direction 2
step=-step;
x_L1_d1 = L1(x_LPC,sysP,step);
%x_dF=x_L1_d1(:,end);%Turning point
xb=x_L1_d1;
xb=[xb(:,1) xb(:,3:end)];

xall=[fliplr(xb) xt];

pl1=true;
    
else 
    pl1=false; %level 1 to run plot script if true
    disp('Level 1 was not computed as per the user request')
end



%The path can be extended to higher dimension level-k by tracking other
%k-extremum points.
%% ----------------------Level 2------------------------------------
prompt="do you want to conduct a level 2 continuation? ";
q=input(prompt,'s');
if q=='y'
%Level 2--Continuation in gamma
fprintf('\n This is the continuation for level 2 \n \n')
fprintf(' ')
prompt="Choose a positive or negative step value: ";
%x_dF=[0;-0.561280292119213;-0.507336003706223;0.00223248564772245;-0.00393174076900323;0;3.28878695322770;0.690984779577617;0.00718200197951464;0.00983498614613223;0;0.179842427559613;0.315710443094361;-0.000159654416876945;0.000145723353895697;-2.88686759750508e-51;-0.126336153545788;0.110179359782989;-0.00254899209413589;-0.00184456579349096;-2.88415210079727e-51;0.586596772841117;-0.785148360082283;0.0100514258510781;-0.00273115448649544;-2.79050794639335e-53;0.105597127321370;-0.00880999910508717;6.76097059658332e-05;0.000172072815521394;1.24441235195691;0.119676714730490];
x_dF=[0;-1.21240194292586;-0.0120529666739570;-0.0267013857110586;-0.000697995361722654;0;1.50778283564525;-0.0893326205881372;0.00677508390313730;-0.000378624638789305;0;0.0493720810034416;0.0914794064024809;1.41481810110253e-05;9.71265276706666e-05;1.77959564569061e-245;-0.0112168620592589;0.624349020090383;-0.00143938221693758;0.0412452165517971;1.76282878463382e-245;-0.0364849634269934;-0.777207247820260;-0.000477031543973740;-0.0104750888829160;-7.36261921111457e-247;0.0474320950645116;-0.0248528771745512;0.000150308270368520;-2.03134733198737e-05;1.49411777562007;0.0230630905743418];
step=input(prompt);
x_L2_d2 = L2(x_dF,sysP,step);

pl2=true; %level 2 variable to run plot script if true

else 
    pl2=false; %level 2 variable to run plot if true
    disp('Level 2 was not computed as per the user request')
end


%----------------------------------------------------
%% Plots
%----------------------MRC---------------------------
x=x_mrc;
figure(1)
subplot(1,2,1)
plot(x(end,:),(sqrt(x(1,:).^2+x(2,:).^2+x(3,:).^2+x(4,:).^2+x(5,:).^2)),'b-','linewidth',1)
xlim([0.5 3.5])
xlabel('$\tilde{F}$','interpreter','latex')
ylabel('$\xi_s$','interpreter','latex')
set(gca,'FontName','times new roman','FontSize',14)
subplot(1,2,2)
plot(x(end,:),(sqrt(x(6,:).^2+x(7,:).^2+x(8,:).^2+x(9,:).^2+x(10,:).^2)),'b-','linewidth',1)
xlim([0.5 3])
xlabel('$\tilde{F}$','interpreter','latex')
ylabel('$\xi_s$','interpreter','latex')
set(gca,'FontName','times new roman','FontSize',14)
%%
%----------------------F-x1-Bifurcation----------------
if pl0
    figure(2)
    plot(x1(end,:),sqrt(x1(1,:).^2+x1(2,:).^2+x1(3,:).^2+x1(4,:).^2+x1(5,:).^2),'Color',[0.4940 0.1840 0.5560],'linewidth',1)
    hold on
    plot(x2(end,:),sqrt(x2(1,:).^2+x2(2,:).^2+x2(3,:).^2+x2(4,:).^2+x2(5,:).^2),'Color',[0.4940 0.1840 0.5560],'linewidth',1)
    plot(x3(end,:),sqrt(x3(1,:).^2+x3(2,:).^2+x3(3,:).^2+x3(4,:).^2+x3(5,:).^2),'Color',[0.4940 0.1840 0.5560],'linewidth',1)
    set(gca,'FontName','times new roman','FontSize',14)
    xlabel('$\tilde{F}$','interpreter','latex')
    ylabel('$\xi_s$','interpreter','latex')
    xlim([0.05,1])
else
    fprintf("\n No variable to plot figure \n ")
end
%%
%-----------------------L1----------------------------
if pl1
    x=xall;
    figure(3)
    plot3(x(end,:),x(end-1,:),(sqrt(x(1,:).^2+x(2,:).^2+x(3,:).^2+x(4,:).^2+x(5,:).^2)),'b-','linewidth',1)
    view(3)
    xlim([0 1])
    xlabel('$\tilde{F}$','interpreter','latex')
    ylabel('$\omega$','interpreter','latex')
    zlabel('$\xi_s$','interpreter','latex')
    set(gca,'FontName','times new roman','FontSize',14)
else
    fprintf("\n No variable to plot figure \n ")
end

%
%%
%-----------------------L2----------------------------
if pl2
    x=x_L2_d2;
    figure(3)
    plot3(x(end,:),x(end-1,:),(sqrt(x(1,:).^2+x(2,:).^2+x(3,:).^2+x(4,:).^2+x(5,:).^2)),'b-','linewidth',1)
    view(2)
    xlim([0 1])
    xlabel('$\tilde{F}$','interpreter','latex')
    ylabel('$\gamma$','interpreter','latex')
    zlabel('$\xi_s$','interpreter','latex')
    set(gca,'FontName','times new roman','FontSize',14)
else
    fprintf("\n No variable to plot figure \n ")
end



%% Optimal parameters determination
% This parameters correspond to the classical arrangement of a nonlinear
% absorber without an inerter. (Nonlinear Den-Hartogs method for equal
% peaks)
    %Mass (kg)
    ms=1;%8.46*10^7;
    eps=0.02;
    mt=eps*ms;
    b=0.02;
%     MassP=[ms,mt,b];
%     
    %Spring coefs (N/m)
    ks=1;%0.85;%7.12*10^7;
    ksn=0.5;
    kt=8*eps*ks*(16+23*eps+9*eps^2+2*(2+eps)*sqrt(4+3*eps))/...\
        (3*(1+eps)^2*(64+80*eps+27*eps^2));
    ktn=2*eps^2*ksn/(1+4*eps);%0.01;
%     SpringP=[ks,ksn,kt,ktn];
%     
    %damping coefs (N-s/m)
    cs=0.002;
    ct=sqrt((kt*mt*(8+9*eps-4*sqrt(4+3*eps)))/(4*(1+eps)));
%%
%The following parameters are generated following the optimal results of
%Joubaneh et al. 
%Nondimensional parameters Joubaneh et. al,   
    ft=1.455;
    fe=1.008;
    uk=0.087;
    Xe=0.188;

% Mass (kg)    
    ms=1;%8.46*10^7;
    mt=0.02;%6.60*10^5;%less then 5% of primary mass 
    b=0.02;%3.30*10^5;
% Spring (N/m)
    ks=1;%0.85;%7.12*10^7;
    %ksn=.0;
    %ktn=0;%7.4074e-04;%0.01;
% Damping  (Ns/m)  
    cs=0.01;
    %ct=0.0;%128;
    %circuit constants

    C=1.02;
    
    wn=sqrt(ks/ms);
    kt=ft^2*wn^2*mt;
    we=fe*wn;
    L=1/(C*we^2);
    kf=sqrt(uk*(L*kt));
    R=Xe*(2*L*we);
    
    kv=kf;%.150;%150;
    
%%
%     wp=sqrt(ks/ms);
%     F0=1;
%     alpha=ksn*F0^2/ks^3
%     zeta=cs/sqrt(ks*ms)
%     beta=kt/ks
%     m=(mt+b)/ms
%     ga=ktn*F0^2/ks^3
%     ld=kf/(F0)
%     sig=R/(L*wp)
%     kap=1/(L*C*wp^2)
%     rho=kv*F0/(ks*L*wp)

