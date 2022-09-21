function [x] = FRF(sysP)
%%finding the periodic solution using the method of harmonic balance

global mu Fw 

y=zeros(15,1);
tic
for kk=1:1
    tic
    %Fw;
    %branch=0;
    if Fw<0.2
        branch=500;
    elseif Fw<0.3
        branch=1000;
    elseif Fw<0.4
        branch=2000;
    elseif Fw<0.5
        branch=3000;
    elseif Fw<0.7
        branch=4000;
    elseif Fw<0.8
        branch=4000;
    else
        branch=10000;
    end
    last_omega=0;
    while last_omega<10
        branch=branch+500;
    mu0=0.5;
    mu1=0.55;
    mu=mu0;
    x0=newton('nondim_temp',y,sysP);
    mu=mu1;
    x1=newton('nondim_temp',x0,sysP);
    [x,p]=branch_follow('nondim_temp',branch,mu0,mu1,x0,x1,sysP);
    last_omega=x(end,end)
    if p=='n'
        break
    end
    end
    %keyboard
kk
toc
end
toc


end

