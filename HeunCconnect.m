% confluent Heun function, a solution of the equation
% HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
% by analytic continuation from point z0, where HeunC(z0) = H0, HeunC'(z0) = dH0,
% to another point z, along the line [z0,z], using a consequence of power expansions
%
% Assumptions:
% z0, z are not 0 or 1
% Im(z0)*Im(z) > 0
%
% Usage:
% [val,dval,err,numb,wrnmsg,R] = HeunCconnect(q,alpha,gamma,delta,epsilon,z,z0,H0,dH0,R0)
%
% R0 is an optional parameter, step size's guess
%
% Returned parameters:
% val is the value of the Heun function
% dval is the value of z-derivative of the Heun function
% err is the estimated error
% numb is the total number of power series terms needed for the evaluation
% wrnmsg is empty if computations are ok
%   otherwise it is a diagnostic message and the function returns val, dval = NaN
% R is the size of the last used step
%
% Oleg V. Motygin, copyright 2018, license: GNU GPL v3
%
% 09 January 2018
%
function [val,dval,err,numb,wrnmsg,R] = HeunCconnect(q,alpha,gamma,delta,epsilon,z,z0,H0,dH0,varargin)

  wrnmsg = '';

  if (z==0)||(z==1)||(z0==0)||(z0==1)||(imag(z)*imag(z0)<0)

    wrnmsg = 'HeunCconnect: assumed that z, z0 are not equal to 0, 1, and Im(z)*Im(z0)>0; ';
    val = NaN; dval = NaN; err = NaN; numb = NaN;

  else

    global Heun_cont_coef Heun_optserterms;
  
    if isempty(Heun_cont_coef) || isempty(Heun_optserterms)
      HeunOpts();
    end
    
    theta = angle(z-z0);
    insearch = true;
    failure = false;
    
    if length(varargin)>0
      R = min(varargin{1},min(abs([z0,z0-1]))*Heun_cont_coef);
    else
      Rmax = 12/(1+abs(epsilon));
      R = min(Rmax,min(abs([z0,z0-1])))*Heun_cont_coef;
    end

    Rtuned = false;
    
    iter = 1;
    
    while ~Rtuned

      if abs(z-z0) <= R
        z1 = z;
      else
        z1 = z0 + R * exp(1i*theta);
      end
      
      [H1,dH1,err,numb,wrnmsg] = HeunCfromZ0(q,alpha,gamma,delta,epsilon,z1,z0,H0,dH0);

      if ~isempty(wrnmsg)
        wrnmsg = strcat('HeunCconnect: problem invoking HeunCfromZ0(',num2str(q),',',num2str(alpha),',', ...
            num2str(gamma),',',num2str(delta),',',num2str(epsilon),',',num2str(z1),',',...
            num2str(z0),',',num2str(H0),',',num2str(dH0),'); warning: ',wrnmsg,'; ');
        failure = true;
        break;      
      end
 
      Rtuned = (err < 5*eps) && (numb < Heun_optserterms) || (iter>5) || (numb<=8);
      
      if ~Rtuned
        R = R / max(err/(5*eps), numb/Heun_optserterms);
      end

      insearch = ~(Rtuned && (z==z1));
      iter = iter+1;
      
    end
    
    z0 = z1; errsum = err; numbsum = numb; H0 = H1; dH0 = dH1;
    
    while insearch && ~failure

      R = min(R,min(abs([z0,z0-1]))*Heun_cont_coef);

      if abs(z-z0) <= R
        z1 = z; insearch = false;
      else
        z1 = z0 + R * exp(1i*theta);
      end
      
      [H0,dH0,err,numb,wrnmsg] = HeunCfromZ0(q,alpha,gamma,delta,epsilon,z1,z0,H0,dH0);

      if ~isempty(wrnmsg)
        wrnmsg = strcat('HeunCconnect: problem invoking HeunCfromZ0(',num2str(q),',',num2str(alpha),',', ...
            num2str(gamma),',',num2str(delta),',',num2str(epsilon),',',num2str(z1),',',...
            num2str(z0),',',num2str(H0),',',num2str(dH0),'); warning: ',wrnmsg,'; ');
        failure = true;
        break;
      end
           
      errsum = errsum + err;
      numbsum = numbsum + numb;

      if insearch
        R = Heun_optserterms * R / (numb + eps);
      end
      
      z0 = z1;
        
    end
  
    numb = numbsum; 
  
    if failure
    
      val = NaN; dval = NaN; err = NaN;
    
    else

      val = H0; dval = dH0; err = errsum;

    end
    
  end
  
end

