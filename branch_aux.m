function y=branch_aux(x,sysP)

global mu tracking_file_name xc arc

% The parameter mu and the solution vector x1 are extracted from the 
% input vector x.
mu=x(end);
x1=x(1:end-1);

% Next the output y is obtained by evaluating the function in the string
% tracking_file_name as y1 and adding the arclength constraint to it. 

y1=feval(tracking_file_name,x1,sysP);
y=[y1;norm([x1;mu]-xc)-arc];

