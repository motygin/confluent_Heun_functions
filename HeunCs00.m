% confluent Heun function, a solution of the equation
% HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
% the second local solution at z=0
% single-valued function with branch-cuts (-\infty,0), (1,+\infty)
% HeunCs(z) = z^(1-gamma)*h(z), where h(0)=1, 
% h'(0)=(-q+(1-gamma)*(delta-epsilon))/(2-gamma) for gamma not equal to 1, 2
% h'(z)/log(z) -> -q+(1-gamma)*(delta-epsilon) as z\to0 for gamma=2
% and
% HeunCs(z) \sim log(z) - q * z * log(z) + ...  as z\to0 for gamma=1
%
% |z| should not exceed the convergency radius 1
%
% Usage:
% [val,dval,err,numb,wrnmsg] = HeunCs00(q,alpha,gamma,delta,epsilon,z)
%
% Returned parameters:
% val is the value of the confluent Heun function
% dval is the value of z-derivative of the confluent Heun function
% err is the estimated error
% numb is the total number of power series terms needed for the evaluation
% wrnmsg is empty if computations are ok
%   otherwise it is a diagnostic message and the function returns val, dval = NaN
%
% Oleg V. Motygin, copyright 2018, license: GNU GPL v3
%
% 09 January 2018
%
function [val,dval,err,numb,wrnmsg] = HeunCs00(q,alpha,gamma,delta,epsilon,z)

  if (abs(z)>=1)

    wrnmsg = 'HeunCs00: z is out of the convergence radius = 1; ';
    val = NaN; dval = NaN; err = NaN; numb = 0;

  elseif (imag(z)==0)&&(real(z)<0)
    
    wrnmsg = 'HeunCs00: z belongs to a possible branch cut; '; 
    val = NaN; dval = NaN; err = NaN; numb = NaN;

  else

    if (abs(gamma-1)<eps)
    
      if ( z==0 )
        val = Inf; dval = Inf;
        err = NaN; numb = 1;
        wrnmsg = '';
      else
        [val,dval,err,numb,wrnmsg] = HeunCs00gamma1(q,alpha,delta,epsilon,z);
      end
  
    else
    
      [H0w,dH0w,errw,numb,wrnmsg] = HeunC00(q+(gamma-1)*(delta-epsilon),alpha+epsilon*(1-gamma),2-gamma,delta,epsilon,z);
      val = z^(1-gamma)*H0w;
      dval = (1-gamma)*z^(-gamma)*H0w + z^(1-gamma)*dH0w;
      if ( isinf(val) || isinf(dval) )
        err = NaN;
      else
        err = abs(z^(1-gamma))*errw;
      end
      
    end
  
  end
  
end

% confluent Heun function, second local solution at z=0, gamma = 1
%
function [val,dval,err,numb,wrnmsg] = HeunCs00gamma1(q,alpha,delta,epsilon,z)
  
  global Heun_klimit;
  
  if isempty(Heun_klimit)
    HeunOpts();
  end

  wrnmsg = '';

  recur0 = @(k,ckm1,ckm2) ...
    (ckm1*z*(-q+(k-1)*(-epsilon+delta+k-1)) + ckm2*z^2*((k-2)*epsilon+alpha))/k^2;

  recur1 = @(k,ckm1,ckm2,dkm0,dkm1,dkm2) recur0(k,ckm1,ckm2) - ...
      (dkm0*2 + dkm1*z*(epsilon/k-delta/k-2+2/k))/k + dkm2*z^2*epsilon/k^2;
      
  L1 = 0; dL1 = 0; ddL1 = 0; dm1 = 0; dm2 = NaN; ckm0 = NaN; ckm1 = 0; ckm2 = 0; 

  L2 = 1; dL2 = 0; ddL2 = 0; skm2 = 0; skm1 = 1;
    
  dsm1 = 0; dsm2 = NaN; skm0 = NaN;
    
  k = 1;

  while ( (k<=Heun_klimit) && ( (dsm2~=dsm1) || (abs(skm0)>eps) || ...
      (dm2~=dm1) || (abs(ckm0)>eps) ) )

    skm0 = recur0(k,skm1,skm2);
    ckm0 = recur1(k,ckm1,ckm2,skm0,skm1,skm2);

    L1 = L1+ckm0; dL1 = dm1+k*ckm0/z; ddL1 = ddL1+k*(k-1)*ckm0/z^2;
    ckm2 = ckm1; ckm1 = ckm0;
    dm2 = dm1; dm1 = dL1;
    
    L2 = L2+skm0; dL2 = dsm1+k*skm0/z; ddL2 = ddL2+k*(k-1)*skm0/z^2;
    skm2 = skm1; skm1 = skm0;
    dsm2 = dsm1; dsm1 = dL2;

    k = k+1;

  end

  numb = k-1;

  val = L1 + log(z) * L2;
  dval = dL1 + log(z) * dL2 + L2/z;
  ddval = ddL1 - L2/z^2 + 2*dL2/z + log(z) * ddL2;

  if ( isinf(val) || isinf(dval) || isnan(val) || isnan(dval) )

    wrnmsg = 'HeunCs00 at gamma=1: failed convergence of recurrence and summation; '; 
    val = NaN; dval = NaN; err = NaN;
      
  else

    if q-alpha*z~=0
  
      val2 = ( z*(z-1)*ddval+(z-1+delta*z+epsilon*z*(z-1))*dval ) / (q-alpha*z);
    
      val3 = ((dL2*epsilon+ddL2)*z^2*log(z)+(dL2*(-epsilon+delta+1)-ddL2)*z*log(z)-dL2*log(z)+ ...
             (dL1*epsilon+ddL1)*z^2+(dL1*(-epsilon+delta+1)+L2*epsilon-ddL1+2*dL2)*z-L2*epsilon+ ...
             L2*delta-2*dL2-dL1) / (q-alpha*z);
             
      err1 = min(abs(val-val2),abs(val-val3));
        
    else
      
      err1 = Inf;
    
    end
      
    if (abs(q-alpha*z)<0.01)||(err1<eps)
    
      err2 = abs(ckm0)*sqrt(numb) + abs(L1)*eps*numb + ...
             abs(log(z)) * ( abs(skm0)*sqrt(numb) + abs(L2)*eps*numb );           
      err =  err2+min(err1,err2);
    
    else
      
      err = err1;
      
    end
       
  end

end
