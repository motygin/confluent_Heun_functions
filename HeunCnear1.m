% confluent Heun function, a solution of the equation
% HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
%
% computation near z=1, by analytic continuation from the point
%
% assumed that z is not equal to 0
%
% computes both the first local solution (see HeunC00),
% and the second local solution (see HeunCs0)
%
% Usage:
% [val1,dval1,err1,val2,dval2,err2,numb,wrnmsg,val1woexp,dval1woexp,err1woexp,val2woexp,dval2woexp,err2woexp] = HeunCnear1(q,alpha,gamma,delta,epsilon,z)
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
%   otherwise it is a diagnostic message and the function returns val, dval = NaN
%
% val1woexp = val1 * exp(epsilon*z)
% dval1woexp = dval1 * exp(epsilon*z)
% err1woexp = err1 * abs(exp(epsilon*z))
%
% val2woexp = val2 * exp(epsilon*z)
% dval2woexp = dval2 * exp(epsilon*z)
% err2woexp = err2 * abs(exp(epsilon*z))
%
% Oleg V. Motygin, copyright 2017-2018, license: GNU GPL v3
%
% 27 December 2017
%
function [val1,dval1,err1,val2,dval2,err2,numb,wrnmsg,val1woexp,dval1woexp,err1woexp,val2woexp,dval2woexp,err2woexp] = HeunCnear1(q,alpha,gamma,delta,epsilon,z)
  
  val1 = NaN; dval1 = NaN; err1 = NaN; numb = 0;
  val2 = NaN; dval2 = NaN; err2 = NaN;
  val1woexp = NaN; dval1woexp = NaN; err1woexp = NaN;
  val2woexp = NaN; dval2woexp = NaN; err2woexp = NaN;

  [m,errJ01,numbJ01,wrnmsgJ01] = HeunCjoin10(q,alpha,gamma,delta,epsilon);
  
  [val1f,dval1f,err1f,numb1f,wrnmsg1f,val1fwoexp,dval1fwoexp,err1fwoexp] = HeunC1(q,alpha,gamma,delta,epsilon,z);
  [val1s,dval1s,err1s,numb1s,wrnmsg1s,val1swoexp,dval1swoexp,err1swoexp] = HeunCs1(q,alpha,gamma,delta,epsilon,z);

  wrnmsg = strcat(wrnmsgJ01,wrnmsg1f,wrnmsg1s);
  numb = numbJ01 + numb1f + numb1s;
  err = errJ01 + err1f + err1s;
  errwoexp = errJ01 + err1fwoexp + err1swoexp;
      
  if isempty(wrnmsg)
      
    val1 = m(1,1)*val1f + m(1,2) * val1s;
    dval1 = m(1,1)*dval1f + m(1,2) * dval1s;
    err1 = abs(m(1,1))*err1f + abs(m(1,2))*err1s;
    
    val2 = m(2,1)*val1f + m(2,2) * val1s;
    dval2 = m(2,1)*dval1f + m(2,2) * dval1s;
    err2 = abs(m(2,1))*err1f + abs(m(2,2))*err1s;
  
    val1woexp = m(1,1)*val1fwoexp + m(1,2) * val1swoexp;
    dval1woexp = m(1,1)*dval1fwoexp + m(1,2) * dval1swoexp;
    err1woexp = abs(m(1,1))*err1fwoexp + abs(m(1,2))*err1swoexp;
    
    val2woexp = m(2,1)*val1fwoexp + m(2,2) * val1swoexp;
    dval2woexp = m(2,1)*dval1fwoexp + m(2,2) * dval1swoexp;
    err2woexp = abs(m(2,1))*err1fwoexp + abs(m(2,2))*err1swoexp;
    
  end
  
end

function [val,dval,err,numb,wrnmsg,valwoexp,dvalwoexp,errwoexp] = HeunC1(q,alpha,gamma,delta,epsilon,z)
    [val,dval,err,numb,wrnmsg,valwoexp,dvalwoexp,errwoexp] = HeunC0(q-alpha,-alpha,delta,gamma,-epsilon,1-z);
    dval = -dval;
    valwoexp = valwoexp * exp(epsilon);
    dvalwoexp = -dvalwoexp * exp(epsilon);
end
   
function [val,dval,err,numb,wrnmsg,valwoexp,dvalwoexp,errwoexp] = HeunCs1(q,alpha,gamma,delta,epsilon,z)
    [val,dval,err,numb,wrnmsg,valwoexp,dvalwoexp,errwoexp] = HeunCs0(q-alpha,-alpha,delta,gamma,-epsilon,1-z);
    dval = -dval;
    valwoexp = valwoexp * exp(epsilon);
    dvalwoexp = -dvalwoexp * exp(epsilon);
end
