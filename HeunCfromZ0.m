% confluent Heun function, a solution of the equation
% HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
% computed at z by power series about Z0 for the given values H(Z0)=H0, H'(Z0)=dH0 
%
% it is assumed that z, Z0 are not equal to 0, 1 and |z-Z0| < min{|Z0|,|Z0-1|}
%
% Usage:
% [val,dval,err,numb,wrnmsg] = HeunCfromZ0(q,alpha,gamma,delta,epsilon,z,Z0,H0,dH0)
%
% Returned parameters:
% val is the value of the confluent Heun function at point z
% dval is the value of z-derivative of the Heun function at point z
% err is the estimated error
% numb is the number of power series terms needed for the evaluation
% wrnmsg is empty if computations are ok
%   otherwise it is a diagnostic message and the function returns val, dval = NaN
%
% Oleg V. Motygin, copyright 2018, license: GNU GPL v3
%
% 09 January 2018
%
function [val,dval,err,numb,wrnmsg] = HeunCfromZ0(q,alpha,gamma,delta,epsilon,z,Z0,H0,dH0)
  
  global Heun_klimit;
  
  if isempty(Heun_klimit)
    HeunOpts();
  end

  R = min(abs(Z0),abs(Z0-1));
  
  wrnmsg = '';

  if (abs(z-Z0)>=R)

    wrnmsg = strcat('HeunCfromZ0: z is out of the convergence radius = ',num2str(R)); val = NaN; dval = NaN; err = NaN; numb = 0;

  elseif ((abs(z-1)<eps) || (abs(Z0-1)<eps))

    wrnmsg = 'HeunCfromZ0: z or Z0 is too close to the singular points'; 
    val = NaN; dval = NaN; err = NaN; numb = 0;

  elseif (z==Z0)

    val = H0; dval = dH0; 
    err = 0; numb = 0;

  else
  
    zeta = z-Z0;
  
    recur = @(k,ckm1,ckm2,ckm3) ...
     (ckm1*zeta*(k-1)*(epsilon*Z0^2+(gamma-epsilon+delta+2*(k-2))*Z0-gamma-k+2)+ ...
        ckm2*zeta^2*((2*(k-2)*epsilon+alpha)*Z0-q+(k-2)*(gamma-epsilon+delta+k-3))+ ...
        ckm3*zeta^3*((k-3)*epsilon+alpha)) / (Z0*(Z0-1)*(1-k)*k);
        
    ckm3 = H0; ckm2 = dH0*zeta; ckm1 = recur(2,ckm2,ckm3,0);

    val = ckm3 + ckm2 + ckm1; vm1 = val; vm2 = NaN;
    dm2 = dH0; dm1 = dH0 + 2*ckm1/zeta; dval = dm1;
    ddval = 2*ckm1/zeta^2; 
  
    k = 3; ckm0 = 1;
    
    while (k<=Heun_klimit) && ( ( vm2~=vm1 ) || ( dm2~=dm1 ) || (abs(ckm0)>eps) )
    
      ckm0 = recur(k,ckm1,ckm2,ckm3);
      val = val + ckm0; dval = dm1 + k*ckm0/zeta;
      ddval = ddval + k*(k-1)*ckm0/zeta^2;
      ckm3 = ckm2; ckm2 = ckm1; ckm1 = ckm0;
      vm2 = vm1; vm1 = val;
      dm2 = dm1; dm1 = dval;
      k=k+1;
      
    end

    numb = k-1;

    if ( isinf(val) || isinf(dval) || isnan(val) || isnan(dval) )
    
      wrnmsg = 'HeunCfromZ0: failed convergence of recurrence and summation; '; 
      val = NaN; dval = NaN; err = NaN;
    
    else
    
      if q-alpha*z~=0
        val2 = ( z*(z-1)*ddval+(gamma*(z-1)+delta*z+epsilon*z*(z-1))*dval ) / (q-alpha*z);
        err1 = abs(val-val2);
      else
        err1 = Inf;
      end
      
      if abs(q-alpha*z)<0.01
        err2 = abs(ckm0) * sqrt(numb) + abs(val) * eps * numb;
        err =  min(err1,err2);
      else
        err = err1;
      end

    end

  end
  
end
