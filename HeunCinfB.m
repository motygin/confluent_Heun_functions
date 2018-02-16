% confluent Heun function, a solution of the equation
% HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
%
% asymptotic expansion at z=infinity,
% the second solution, including exponential factor
%
% Usage:
% [val,dval,err,numb,wrnmsg,valwoexp,dvalwoexp,errwoexp] = HeunCinfB(q,alpha,gamma,delta,epsilon,z)
%
% Returned parameters:
% val is the value of the Heun function
% dval is the value of z-derivative of the Heun function
% err is the estimated error
% numb is the number of the summed series terms
% wrnmsg is a warning message:
%   it is empty if computations are ok
%   otherwise it is a diagnostic message and the function returns val, dval = NaN
%
% the following values are also returned:
%
% valwoexp = val * exp(epsilon*z)
% dvalwoexp = dval * exp(epsilon*z)
% errwoexp = err * abs(exp(epsilon*z))
%
% Oleg V. Motygin, copyright 2017-2018, license: GNU GPL v3
%
% 20 December 2017
%
function [val,dval,err,numb,wrnmsg,valwoexp,dvalwoexp,errwoexp] = HeunCinfB(q,alpha,gamma,delta,epsilon,z)

  [val0,dval0,err0,numb,wrnmsg] = HeunCinfA(q-epsilon*gamma,alpha-epsilon*(gamma+delta),gamma,delta,-epsilon,z);
  
  valwoexp = val0;
  dvalwoexp = -epsilon * val0 + dval0;
  errwoexp = err0;
  
  val = exp(-epsilon*z) * valwoexp;
  dval = exp(-epsilon*z) * dvalwoexp;
  err = abs(exp(-epsilon*z)) * errwoexp;

end

    
