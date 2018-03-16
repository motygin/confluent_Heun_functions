% confluent Heun function, a solution of the equation
% HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
%
% computation for sufficiently large |z|, by analytic continuation from infinity
%
% computes both the first at z=0 local solution (see HeunC00) and the second at z=0 local solution (see HeunCs0)
%
% Usage:
% [val1,dval1,err1,val2,dval2,err2,numb,wrnmsg,val1woexp,dval1woexp,err1woexp,val2woexp,dval2woexp,err2woexp] = HeunCfaraway(q,alpha,gamma,delta,epsilon,z)
%
% Returned parameters:
% val1 is the value of the Heun function, growing from the first local solution at z=0
% dval1 is the value of z-derivative of the Heun function
% err1 is the estimated error
% val2 is the value of the Heun function, growing from the second local solution at z=0
% dval2 is the value of z-derivative of the Heun function
% err2 is the estimated error
% numb is the number of power series terms needed for the evaluation
% wrnmsg is a warning message:
%   it is empty if computations are ok
%   otherwise it is a diagnostic message and the function returns val*, dval* = NaN
%
% val1woexp = val1 * exp(epsilon*z)
% dval1woexp = dval1 * exp(epsilon*z)
% err1woexp = err1 * abs(exp(epsilon*z))
%
% val2woexp = val2 * exp(epsilon*z)
% dval2woexp = dval2 * exp(epsilon*z)
% err2woexp = err2 * abs(exp(epsilon*z))
%
% Oleg V. Motygin, copyright 2018, license: GNU GPL v3
%
% 16 March 2018
%
function [val1,dval1,err1,val2,dval2,err2,numb,wrnmsg,val1woexp,dval1woexp,err1woexp,val2woexp,dval2woexp,err2woexp] = HeunCfaraway(q,alpha,gamma,delta,epsilon,z)
  
  val1 = NaN; dval1 = NaN; err1 = NaN; numb = 0;
  val2 = NaN; dval2 = NaN; err2 = NaN;
  val1woexp = NaN; dval1woexp = NaN; err1woexp = NaN;
  val2woexp = NaN; dval2woexp = NaN; err2woexp = NaN;

  hs = sign(imag(z));
  
  [C0A,CsA,errJA,numbJA,wrnmsgJA] = HeunCjoin0infA(q,alpha,gamma,delta,epsilon,hs*pi/2);
  [C0B,CsB,errJB,numbJB,wrnmsgJB] = HeunCjoin0infB(q,alpha,gamma,delta,epsilon,hs*pi/2);
  
  [R,N] = findR();
  infpt = 2 * max(1,R/(abs(eps)+abs(epsilon))) * z/abs(z);

  if abs(z)>abs(infpt)
  
    [valA,dvalA,errA,numbA,wrnmsgA] = HeunCinfA(q,alpha,gamma,delta,epsilon,z);
    [valB,dvalB,errB,numbB,wrnmsgB] = HeunCinfA(q-epsilon*gamma,alpha-epsilon*(gamma+delta),gamma,delta,-epsilon,z);
    
    wrnmsg = strcat(wrnmsgJA,wrnmsgJB,wrnmsgA,wrnmsgB);
    numb = numbJA + numbJB + numbA + numbB;
    err = errJA + errJB + errA + errB;
    
  else  
    
    [valinfA,dvalinfA,errinfA,numbinfA,wrnmsginfA] = HeunCinfA(q,alpha,gamma,delta,epsilon,infpt);
    [valA,dvalA,errA,numbA,wrnmsgA] = HeunCconnect(q,alpha,gamma,delta,epsilon,z,infpt,valinfA,dvalinfA);
    
    [valinfB,dvalinfB,errinfB,numbinfB,wrnmsginfB] = HeunCinfA(q-epsilon*gamma,alpha-epsilon*(gamma+delta),gamma,delta,-epsilon,infpt);
    [valB,dvalB,errB,numbB,wrnmsgB] = HeunCconnect(q-epsilon*gamma,alpha-epsilon*(gamma+delta),gamma,delta,-epsilon,z,infpt,valinfB,dvalinfB);
    
    wrnmsg = strcat(wrnmsgJA,wrnmsgJB,wrnmsginfA,wrnmsgA,wrnmsginfB,wrnmsgB);
    numb = numbJA + numbJB + numbinfA + numbA + numbinfB + numbB;
    err = errJA + errJB + errinfA + errA + errinfB + errB;

  end  
      
  if isempty(wrnmsg)
      
    m = inv([C0A,CsA;C0B,CsB]);
    
    co = cond(m);

    val1 = m(1,1) * valA + m(1,2) * exp(-epsilon*z) * valB;
    dval1 = m(1,1) * dvalA + m(1,2) * exp(-epsilon*z) * (-epsilon * valB + dvalB);
    err1 = co * (abs(m(1,1)) * errA + abs(m(1,2)) * abs(exp(-epsilon*z)) * errB);
    
    val2 = m(2,1) * valA + m(2,2) * exp(-epsilon*z) * valB;
    dval2 = m(2,1) * dvalA + m(2,2) * exp(-epsilon*z) * (-epsilon * valB + dvalB);
    err2 = co * (abs(m(2,1)) * errA + abs(m(2,2)) * abs(exp(-epsilon*z)) * errB);
    
    val1woexp = m(1,1) * exp(epsilon*z) * valA + m(1,2) * valB;
    dval1woexp = m(1,1) * exp(epsilon*z) * dvalA + m(1,2) * (-epsilon * valB + dvalB);
    err1woexp = co * (abs(m(1,1)) * abs(exp(epsilon*z)) * errA + abs(m(1,2)) * errB);
    
    val2woexp = m(2,1) * exp(epsilon*z) * valA + m(2,2) * valB;
    dval2woexp = m(2,1) * exp(epsilon*z) * dvalA + m(2,2) * (-epsilon * valB + dvalB);
    err2woexp = co * (abs(m(2,1)) * abs(exp(epsilon*z)) * errA + abs(m(2,2)) * errB);
    
  end
    
end
