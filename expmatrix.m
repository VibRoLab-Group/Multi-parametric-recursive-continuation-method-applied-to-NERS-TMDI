function Ab=expmatrix(t,h,y1,sysP,omega1,F)
%In this function set of linear odes with periodic coeffecients
%have been produced for the calculation of floquet multiplier.

 global alpha zeta

%     %system parameters
%     alpha=1;
%     zeta=0.01;
    
    %tuned parameters
    b=sysP(1);
    mu=sysP(2);
    ld=sysP(3);
    ga=sysP(4);
    sig=sysP(5);
    kap=sysP(6);
    rho=sysP(7);
    th=sysP(8);
    
omega=omega1;
y=y1;
%keyboard
x_s1=y(1)+y(2)*cos(omega*t)+y(3)*sin(omega*t)+y(4)*cos(3*omega*t)+y(5)*sin(3*omega*t);
x_s3=y(6)+y(7)*cos(omega*t)+y(8)*sin(omega*t)+y(9)*cos(3*omega*t)+y(10)*sin(3*omega*t);
%x_s5=y(11)+y(12)*cos(omega*t)+y(13)*sin(omega*t)+y(14)*cos(3*omega*t)+y(15)*sin(3*omega*t);
%keyboard
% f=[0 -(-x(2,1)) 0 0 0 0;...
%     -x(1,1)*(b + 3*alpha*x_s1^2+3*ga*(x_s1 - x_s3)^2 + 1) -x(2,1)*zeta  x(3,1)*(b+3*ga*(x_s1-x_s3)^2) 0 0 0;...
%     0 0 0 -(-x(4,1)) 0 0;...
%     x(1,1)*(b/mu+(3*ga*(x_s1-x_s3)^2)/mu) 0 -x(3,1)*(b/mu+(3*ga*(x_s1-x_s3)^2)/mu) 0 0 -ld*x(6,1)/mu ;...
% 	0 0 0 0 0 -(-x(6,1));...
%     0 0 0 rho*x(4,1) -kap*x(5,1) -sig*x(6,1)];

f=[0 1 0 0 0 0;...
    -(b + 3*alpha*x_s1^2+3*ga*(x_s1 - x_s3)^2 + 1) -zeta-th  (b+3*ga*(x_s1-x_s3)^2) th 0 0;...
    0 0 0 1 0 0;...
    (b/mu+(3*ga*(x_s1-x_s3)^2)/mu) th/mu -(b/mu+(3*ga*(x_s1-x_s3)^2)/mu) -th/mu 0 -ld/mu ;...
	0 0 0 0 0 1;...
    0 0 0 rho -kap -sig];
% f=[0 1 0 0 ;...
%     -(b + 3*alpha*x_s1^2+3*ga*(x_s1 - x_s3)^2 + 1) -zeta  (b+3*ga*(x_s1-x_s3)^2) 0;...
%     0 0 0 1 ;...
%     (b/mu+(3*ga*(x_s1-x_s3)^2)/mu) 0 -(b/mu+(3*ga*(x_s1-x_s3)^2)/mu) 0];
%keyboard

Ab=expm(f*h);
end