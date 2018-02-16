% for confluent Heun function, a solution of the equation
% HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
%
% HeunCjoin0infB finds connection coefficients C0, Cs, such that
% C0 * HeunC00(z) + Cs * HeunCs00(z) analytically continues to the second, 
% including exponential factor solution at infinity \exp(i\theta) \infty
% (see HeunCinfB)
%
% Usage:
% [C0,Cs,err,numb,wrnmsg] = HeunCjoin0infB(q,alpha,gamma,delta,epsilon,theta)
%
% Returned parameters:
% C0, Cs are the connection coefficients
% err is the estimated error
% numb is the number of power series terms needed for the evaluation
% wrnmsg is a warning message:
%   it is empty if computations are ok
%   otherwise it is a diagnostic message and the function returns C0, Cs = NaN
%
% Oleg V. Motygin, copyright 2018, license: GNU GPL v3
%
% 15 February 2018
%
function [C0,Cs,err,numb,wrnmsg] = HeunCjoin0infB(q,alpha,gamma,delta,epsilon,theta)
  
  [C0,Cs,err,numb,wrnmsg] = HeunCjoin0infA(q-epsilon*gamma,alpha-epsilon*(gamma+delta),gamma,delta,-epsilon,theta,'woexp');
  
end
