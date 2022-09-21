function x = L1(x_interest,sysP,step)
%%finding the periodic solution using the method of harmonic balance
clc
%clear all
close all
global mu Fw count

%ga=sysP(4);

x0=x_interest(1:end-1)';
phi0=[1.00527873511279e-05,-0.129311522027891,0.0840882766970096,0.00104000886528385,-0.00480988693351045,6.81803296626275e-06,0.664108575557058,-0.723380592232881,0.00847242449023657,-0.000888433735559738,-2.87194440146902e-05,0.108585484376654,-0.000273171989298249,7.98198223145913e-05,0.000577355643005017];
w0=x_interest(end);

y_uncalibrated=[x0,phi0,w0]';

%determine phi0
count=1;
phi0=nondim_temp2(y_uncalibrated,sysP);

count=2;
y=[x0,phi0,w0]';


tic
for kk=1:1
    tic
    branch=30000;
    mu0=Fw;
    mu1=mu0+step;
    mu=mu0;
    x0=newton('nondim_temp2',y,sysP);
    keyboard
    x1=x0;
    [x,p]=branch_follow2('nondim_temp2',branch,mu0,mu1,x0,x1,sysP);
    x=[x(:,1) x(:,3:end)];
    if p=='n'
        break
    end
kk
toc
end
toc
end

