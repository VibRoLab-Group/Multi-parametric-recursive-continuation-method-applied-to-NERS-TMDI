function [x0,c]=newton(fun,x0,sysP)

% This is the code to obtain roots of nonlinear equations using the
% Newton-Raphson method.

% The maximum number of iterations and the tolerance 
% for locating a root is specified here.  
nmax=100;
epsil=1e-6; % tolerance 

n=0; % The step counter is set initially to zero. 
%keyboard
f0=feval(fun,x0,sysP); % The value of the function at the initial guess x0 is evaluated here. 
norm(f0);
% The derivative matrix corresponding to the desired function is required
% for the newton method. In this code, we generate this matrix numerically.
% For this, we require to specify the size of the matrix and the time-step
% to be used for numerical determination of the derivative matrix.

m=length(x0); % size of the derivative matrix.
h=1e-5;       % time-step for evaluating the derivative matrix. 
E=eye(m);     % E is the identity matrix useful for specifying initial conditions 
              % for the numerical evaluation of the derivative matrix.
c=1;          % c is the handle determining of the roots have converged or not.

% We check if the function values are less than the desired tolerance or
% not. This is the check of the convergence of the roots.
%keyboard
while (norm(f0)>epsil*max([1,norm(x0)]))*(n<nmax)
   n=n+1;  % The step counter is incrementd by one.
   D=E;    % The initial derivative matrix is set to the identity matrix.
   % In the next step the derivative matrix is obtained numerically.
   for k=1:m
      D(:,k)=feval(fun,x0+h*E(:,k),sysP)-f0;
   end
   D=D/h;

    % Calculate
   
   % The next guess for the root is evaluated below and the function is
   % evaluated at the new guess.
   x0=x0-D\f0*0.5
   f0=feval(fun,x0,sysP);
end
% Next we check if the roots have converged or not by checking the value of
% the step counter n. 
if n==nmax
    c=0;
   x0=inf*x0;
   disp('did not converge')
end
