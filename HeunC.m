% confluent Heun function, the first local solution of the equation
% HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
% at z=0 such that HeunC(0)=1 and
% HeunC'(0)=-q/gamma when gamma is not equal to 0
% HeunC'(z)/log(z) -> -q as z->0 when gamma = 0
%
% computed by a consequence of power expansions with improvements near points z=1 and z=\infty
%
% it is assumed that z does not belong to the branch-cut [1,\infty)
%
% Usage:
% [val,dval,err,numb,wrnmsg,valwoexp,dvalwoexp,errwoexp] = HeunC(q,alpha,gamma,delta,epsilon,z)
%
% Returned parameters:
% val is the value of the Heun function
% dval is the value of z-derivative of the Heun function
% err is the estimated error
% numb is the number of power series terms needed for the evaluation
% wrnmsg is a warning message:
%   it is empty if computations are ok
%   otherwise it is a diagnostic message and the function returns val, dval = NaN
%
% valwoexp = val * exp(epsilon*z)
% dvalwoexp = dval * exp(epsilon*z)
% errwoexp = err * abs(exp(epsilon*z))
%
% Oleg V. Motygin, copyright 2018, license: GNU GPL v3
%
% 26 January 2018
%
function [val,dval,err,numb,wrnmsg,valwoexp,dvalwoexp,errwoexp] = HeunC(q,alpha,gamma,delta,epsilon,z)
  
  global Heun_proxco Heun_proxcoinf_rel;
  
  if isempty(Heun_proxco) || isempty(Heun_proxcoinf_rel)
    HeunOpts();
  end

  if (imag(z)==0)&&(real(z)>=1)
    wrnmsg = 'HeunC: z belongs to the branch-cut [1,\infty); ';
    val = NaN; dval = NaN; err = NaN; numb = 0;
    valwoexp = NaN; dvalwoexp = NaN; errwoexp = NaN;
  
  else
    
    [R,N] = findR();
    
    if (abs(z-1)<Heun_proxco)

      [val,dval,err,val2,dval2,err2,numb,wrnmsg,valwoexp,dvalwoexp,errwoexp,val2woexp,dval2woexp,err2woexp] = HeunCnear1(q,alpha,gamma,delta,epsilon,z);

    elseif (abs(epsilon)>1/2)&&(abs(q)<2.5)&&(abs(z)>Heun_proxcoinf_rel*R/(abs(eps)+abs(epsilon)))
  
      [val,dval,err,val2,dval2,err2,numb,wrnmsg,valwoexp,dvalwoexp,errwoexp,val2woexp,dval2woexp,err2woexp] = HeunCfaraway(q,alpha,gamma,delta,epsilon,z);
  
    else
      
      [val,dval,err,numb,wrnmsg,valwoexp,dvalwoexp,errwoexp] = HeunC0(q,alpha,gamma,delta,epsilon,z);
    
    end
    
  end

end
