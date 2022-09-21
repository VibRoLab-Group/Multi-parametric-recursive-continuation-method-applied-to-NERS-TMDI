function [x,p]=branch_follow(fname,nsteps,mu0,mu1,x0,x1,sysP)
% The aim of this program is to follow solution branches to systems of nonlinear
% equations with one free parameter. 

% You start with a file fname. This file uses a global parameter (called mu).
% Once mu is defined, given x, this file returns f(x). Note: x is n-dimensional.
% In the space of x, there  is a one-parameter family (or curve) of solutions
% parameterized by mu. Two nearby points on this curve (mu0,x0) and (mu1,x1) are 
% initially specified (found manually with trial and error). Our aim is to generate
% a sequence of points (mu2,x2), (mu3,x3) etc. The distance between each point
% and the next one is to be kept fixed as we move along the curve. Some auxiliary
% files will be needed, and called as appropriate.

global mu tracking_file_name xc arc Fw

% We specify the file name containing the function whose roots are to be
% tracked. This would be used for augmenting the system equations with the 
% arclength constraint in the file branch_aux.
tracking_file_name=fname;

% In the next portion, we extend the solution vector x by 1 by adding
% the varibale mu to it. This will be useful in fixing the arclength of the
% curve to be tracked.

x0=[x0;mu0]; % The first point on the curve.
xc=[x1;mu1]; % The second point on the curve.

% Below we initialize the solution vector x which will be storing the
% extended solutions for the problem.
x=[x0,xc];

% The arclength is fixed based on the first two initial conditions that are
% specifed as inputs to this file.
arc=norm(x0-xc);

% k is the counter number for continuation and c is the handle which
% specifies if the newton method used for getting the roots of the
% equations has converged or not. If the method does not converge, we will
% have c=0. We initialize it to be c=1.
k=1;
c=1;
track=false;
p='';

while (k<nsteps)*c
   xg=2*xc-x0;  % extrapolation to get the initial guess for the next point.
   [xx,c]=newton('branch_aux',xg,sysP); % The roots are obtained using the newton method here.
   if c
      k=k+1; % If we get a converged root we increment the continuation step counter here.
      x0=xc; % The second point now becomes the first point for the next step.
      
      xc=xx; % The new root is defined as the second point here.
      
      x=[x,xx]; % The new root is added to the vector x for storage.
      save result x;
      
      %check stability
        y1=xx(1:10);
        n1=200;
        omega1=xx(end);
        tf=2*pi/(omega1);
        h=tf/n1;
        flmu=temp(h,n1,y1,sysP,omega1,Fw);
        eflmu=eig(flmu);
        eflmu(eflmu<0)=-eflmu(eflmu<0);
        if (max(eflmu)<1.0) && (track==false)
            floquet_multiplier=max(eig(flmu));
            track=false;
        elseif (max(eflmu)>1.0) && (track==false)
            track=true;
            xx(end)
            floquet_multiplier=max(eig(flmu))
            if isreal(floquet_multiplier) && (floquet_multiplier>0) 
                disp("LPC")
                
                prompt="continue? type 'y' for yes or 'n' for no: ";
                p=input(prompt,'s')
                if p=='n'
                    disp('see you next time')
                    break
                end   
            else 
                disp("NS")
                                prompt="continue? type 'y' for yes or 'n' for no: ";
                p=input(prompt,'s')
                if p=='n'
                    disp('see you next time')
                    break
                end   
            end
        elseif (max(eflmu)<1.0) && (track==true)
            floquet_multiplier=max(eig(flmu))
            track=false; 
            prompt="continue? type 'y' for yes or 'n' for no: ";
            p=input(prompt,'s')
            if p=='n'
            disp('see next time')
            break
            end   
        end


      % The below step is just for display and can be turned off as well.
      % xx(2) % The first componenet of the root is displayed here so that we can find if the
            % desired trajectory is being tracked.
   end
end

   