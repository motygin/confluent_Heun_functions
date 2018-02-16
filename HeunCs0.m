% confluent Heun function, a solution of the equation
% HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
% the second local solution at z=0 (see HeunCs00)
%
% computed by a consequence of power expansions
%
% it is assumed that z does not belong to the branch-cuts (-\infty,0], [1,\infty)
%
% Usage:
% [val,dval,err,numb,wrnmsg,valwoexp,dvalwoexp,errwoexp] = HeunCs0(q,alpha,gamma,delta,epsilon,z)
%
% Returned parameters:
% val is the value of the Heun function
% dval is the value of z-derivative of the Heun function
% err is the estimated error
% numb is the number of power series terms needed for the evaluation
% wrnmsg is a warning message:
%   it is empty if computations are ok
%   otherwise it is a diagnostic message and the function returns val*, dval* = NaN
%
% valwoexp = val * exp(epsilon*z)
% dvalwoexp = dval * exp(epsilon*z)
% errwoexp = err * abs(exp(epsilon*z))
%
% Oleg V. Motygin, copyright 2018, license: GNU GPL v3
%
% 15 February 2018
%
function [val,dval,err,numb,wrnmsg,valwoexp,dvalwoexp,errwoexp] = HeunCs0(q,alpha,gamma,delta,epsilon,z)
  
  global Heun_cont_coef;
  
  if isempty(Heun_cont_coef)
    HeunOpts();
  end
  
  wrnmsg = '';

  if (imag(z)==0)&&((real(z)>=1)||(real(z)<=0))
    wrnmsg = 'HeunCs0: z belongs to the possible branch-cuts (-\infty,0], [1,\infty); ';
    val = NaN; dval = NaN; err = NaN; numb = 0;
    valwoexp = NaN; dvalwoexp = NaN; errwoexp = NaN;

  else
  
    expgrow = real(-epsilon*z)>0;

    if expgrow
    
      q = q - epsilon * gamma;
      alpha = alpha - epsilon * (gamma+delta);
      epsilon = -epsilon;

    end  
    
    if (abs(z)<Heun_cont_coef)

      [val,dval,err,numb,wrnmsg] = HeunCs00(q,alpha,gamma,delta,epsilon,z);

    elseif (real(z)>1)&&((imag(z)>0)&&(imag(z)<=real(z))||(imag(z)<0)&&(imag(z)>=-real(z)))
  
      z1 = 1 + sign(imag(z))*1j; z0 = Heun_cont_coef * z1;
      [H0,dH0,err0,numb0,wrnmsg0] = HeunCs00(q,alpha,gamma,delta,epsilon,z0);
      [H1,dH1,err1,numb1,wrnmsg1,R] = HeunCconnect(q,alpha,gamma,delta,epsilon,z1,z0,H0,dH0);
      [val,dval,err2,numb2,wrnmsg2] = HeunCconnect(q,alpha,gamma,delta,epsilon,z,z1,H1,dH1,R);
      numb = numb0 + numb1 + numb2;
      err = err0 + err1 + err2;
      wrnmsg = strcat(wrnmsg0,wrnmsg1,wrnmsg2);
  
    else
      
      z0 = Heun_cont_coef*z/abs(z);
      [H0,dH0,err0,numb0,wrnmsg0] = HeunCs00(q,alpha,gamma,delta,epsilon,z0);
      [val,dval,err1,numb1,wrnmsg1] = HeunCconnect(q,alpha,gamma,delta,epsilon,z,z0,H0,dH0);
      numb = numb0 + numb1;
      err = err0 + err1;
      wrnmsg = strcat(wrnmsg0,wrnmsg1);

    end
    
    if expgrow
    
      valwoexp = val; 
      dvalwoexp = epsilon*val + dval; 
      errwoexp = err;
      val = valwoexp * exp(epsilon*z);
      dval = dvalwoexp * exp(epsilon*z);
      err = err * abs(exp(epsilon*z));

    else
     
      valwoexp = val * exp(epsilon*z);
      dvalwoexp = dval * exp(epsilon*z);
      errwoexp = err * abs(exp(epsilon*z));
      
    end
    
  end

end

