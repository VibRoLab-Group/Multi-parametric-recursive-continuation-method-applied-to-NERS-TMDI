% function flmu = temp(h,n1,N,modulation,omegam,kappa,xi,beta,psi,n,omega0,v)
% %(h,n1,N,modulation,omegam,kappa,xi,beta,psi,n,omega0,v)
% flmu=eye(29);
%        for k=1:n1
%            %keyboard
%             t=(k-1)*h;
%             temp=expmatrix(t,h,N,modulation,omegam,kappa,xi,beta,psi,n,omega0,v);
%             %the aove function evalutes exponetial of a matrix
%             flmu=flmu*temp;
%             %keyboard
%        end
% 
% end


function flmu = temp(h,n1,y1,sysP,omega1,F)
%
flmu=eye(6);
       for k=1:n1
           %keyboard
            t=(k-1)*h;
            temp=expmatrix(t,h,y1,sysP,omega1,F);
            %the above function evalutes exponential of a matrix
            flmu=flmu*temp;
            %keyboard
       end

end
