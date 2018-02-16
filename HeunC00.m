% confluent Heun function, the first local solution of the equation
% HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
% at z=0 such that HeunC(0)=1 and
% HeunC'(0)=-q/gamma when gamma is not equal to 0
% HeunC'(z)/log(z) -> -q as z->0 when gamma = 0
%
% |z| should not exceed the convergency radius 1
%
% Usage:
% [val,dval,err,numb,wrnmsg] = HeunC00(q,alpha,gamma,delta,epsilon,z,aux)
%
% aux is an optional parameter, it only works for gamma = 0, -1, -2, ...
% if aux = "yes" then the function is computed as a combination with the second
% solution so as to satisfy
% HeunC00(q,alpha,gamma,delta,varepsilon;z)=
% exp(-\varepsilon z)[HeunC00(q-epsilon*gamma,alpha-varepsilon(gamma+delta),gamma,delta,-varepsilon;z)+
% + A * HeunCs00(q-epsilon*gamma,alpha-varepsilon(gamma+delta),gamma,delta,-varepsilon;z)]
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
% Oleg V. Motygin, copyright 2018, license: GNU GPL v3
%
% 15 February 2018
%
function [val,dval,err,numb,wrnmsg] = HeunC00(q,alpha,gamma,delta,epsilon,z,varargin)
  
  if (abs(z)>=1)

    wrnmsg = 'HeunC00: z is out of the convergence radius = 1; ';
    val = NaN; dval = NaN; err = NaN; numb = NaN;

  else
  
    isnonpositiveintegergamma = abs(ceil(gamma-5*eps)+abs(gamma))<5*eps;

    if isnonpositiveintegergamma
    
      [val,dval,err,numb,wrnmsg] = HeunC00log(q,alpha,gamma,delta,epsilon,z);
      
      if length(varargin)>0
        if strcmp(varargin{1},'woexp')
      
          co = findcoef4HeunCs(q,alpha,gamma,delta,epsilon);
          [vals,dvals,errs,numbs,wrnmsgs] = HeunCs00(q,alpha,gamma,delta,epsilon,z);

          val = val - co * vals;
          dval = dval - co * dvals;
          err = err + abs(co) * errs;
          numb = numb + numbs;
          wrnmsg = strcat(wrnmsg,wrnmsgs);
        end
      end

    else
    
      [val,dval,err,numb,wrnmsg] = HeunC00gen(q,alpha,gamma,delta,epsilon,z);

    end

  end

end


% confluent Heun function, gamma is not equal to 0, -1, -2, ...
%
function [val,dval,err,numb,wrnmsg] = HeunC00gen(q,alpha,gamma,delta,epsilon,z)
  
  global Heun_klimit;
  
  if isempty(Heun_klimit)
    HeunOpts();
  end

  wrnmsg = '';
  
  if (z==0)
  
    val = 1; dval = -q/gamma;
    err = 0; numb = 1;
  
  else

    recur = @(k,ckm1,ckm2) ...
      (ckm1*z*(-q+(k-1)*(gamma-epsilon+delta+k-2)) + ...
       ckm2*z^2*((k-2)*epsilon+alpha))/(k*(gamma+k-1));

    ckm2 = 1; ckm1 = -z*q/gamma;

    val = ckm2+ckm1; vm1 = val; vm2 = NaN;
    dval = -q/gamma; dm1 = dval; dm2 = NaN;
    ddval = 0;
    
    k = 2; ckm0 = 1;

    while ( (k<=Heun_klimit) && ( (vm2~=vm1) || (dm2~=dm1) || (abs(ckm0)>eps) ) )
      
      ckm0 = recur(k,ckm1,ckm2);
      val = val + ckm0; dval = dm1 + k*ckm0/z;
      ddval = ddval + k*(k-1)*ckm0/z^2;
      ckm2 = ckm1; ckm1 = ckm0;
      vm2 = vm1; vm1 = val;
      dm2 = dm1; dm1 = dval;
      k = k+1;
    
    end

    numb = k-1;

    if ( isinf(val) || isinf(dval) || isnan(val) || isnan(dval) )
    
      wrnmsg = 'HeunC00: failed convergence of recurrence and summation; '; 
      val = NaN; dval = NaN; err = NaN;
      
    else
    
      if q-alpha*z~=0
        val2 = ( z*(z-1)*ddval+(gamma*(z-1)+delta*z+epsilon*z*(z-1))*dval ) / (q-alpha*z);
        err1 = abs(val-val2);
      else
        err1 = Inf;
      end
      
      if abs(q-alpha*z)<0.01
        err2 = abs(ckm0) * sqrt(numb) + eps * numb * abs(val);
        err =  min(err1,err2);
      else
        err = err1;
      end
      
    end

  end

end

% confluent Heun function, gamma = 0, -1, -2, ...
%
function [val,dval,err,numb,wrnmsg] = HeunC00log(q,alpha,gamma,delta,epsilon,z)
  
  global Heun_klimit;
  
  if isempty(Heun_klimit)
    HeunOpts();
  end
  
  wrnmsg = '';

  if (z==0)
  
    val = 1; 
    err = 0; numb = 1;
  
    if (abs(gamma)<eps)
      dval = Inf;
    else
      dval = -q/gamma;
    end

  else
  
    N = round(1-gamma);
  
    recur0 = @(k,ckm1,ckm2) ...
      (ckm1*z*(-q+(k-1)*(gamma-epsilon+delta+k-2)) + ...
       ckm2*z^2*((k-2)*epsilon+alpha))/(k*(gamma+k-1));

    recur1 = @(k,ckm1,ckm2,dkm0,dkm1,dkm2) recur0(k,ckm1,ckm2) + ...
      (-dkm0*(gamma+2*k-1)+dkm1*z*(gamma-epsilon+delta+2*k-3)+dkm2*z^2*epsilon)/(k*(gamma+k-1));   

    L1 = 1; dL1 = 0; ddL1 = 0; ckm0 = 1; ckm1 = 1; ckm2 = 0;
    
    for k=1:N-1
      ckm0 = recur0(k,ckm1,ckm2);
      
      L1 = L1+ckm0; dL1 = dL1+k*ckm0/z; ddL1 = ddL1+k*(k-1)*ckm0/z^2;
      ckm2 = ckm1; ckm1 = ckm0; 
    end
       
    sN = (ckm1*z*(q+gamma*(delta-epsilon-1)) + ...
         ckm2*z^2*(epsilon*(gamma+1)-alpha))/(gamma-1);

    L2 = 0; dL2 = 0; ddL2 = 0; dm1 = dL2; dm2 = NaN;
    ckm1 = 0; ckm2 = ckm0; 

    L3 = sN; skm2 = 0; skm1 = sN;
    dL3 = N*sN/z; ddL3 = N*(N-1)*sN/z^2; 
    dsm1 = dL3; dsm2 = NaN; skm0 = NaN;
    
    k = N+1;

    while ( (k<=Heun_klimit) && ( (dsm2~=dsm1) || (abs(skm0)>eps) || ...
      (dm2~=dm1) || (abs(ckm0)>eps) ) )

      skm0 = recur0(k,skm1,skm2);
      ckm0 = recur1(k,ckm1,ckm2,skm0,skm1,skm2);

      L2 = L2+ckm0; dL2 = dm1+k*ckm0/z; ddL2 = ddL2+k*(k-1)*ckm0/z^2;

      ckm2 = ckm1; ckm1 = ckm0;
      dm2 = dm1; dm1 = dL2;

      L3 = L3+skm0; dL3 = dsm1+k*skm0/z; ddL3 = ddL3+k*(k-1)*skm0/z^2;
      skm2 = skm1; skm1 = skm0;
      dsm2 = dsm1; dsm1 = dL3;

      k = k+1;

    end
    
    numb = k-1;

    val = L1 + L2 + log(z) * L3;
    dval = dL1 + dL2 + log(z) * dL3 + L3/z;
    ddval = ddL1 + ddL2 - L3/z^2 + 2*dL3/z + log(z) * ddL3;

    if ( isinf(val) || isinf(dval) || isnan(val) || isnan(dval) )

      wrnmsg = 'HeunL00log: failed convergence of recurrence and summation; '; 
      val = NaN; dval = NaN; err = NaN;
      
    else
    
      if q-alpha*z~=0
    
        val2 = ( z*(z-1)*ddval+(gamma*(z-1)+delta*z+epsilon*z*(z-1))*dval ) / (q-alpha*z);

        val3 = ((dL3*epsilon+ddL3)*z^2*log(z)+(dL3*(gamma-epsilon+delta)-ddL3)*z*log(z)-dL3*gamma*log(z)+ ...
          (epsilon*(dL2+dL1)+ddL2+ddL1)*z^2+((dL1+dL2)*(gamma-epsilon+delta)+L3*epsilon-ddL2-ddL1+2*dL3)*z+...
          L3*(1-gamma)/z-(dL1+dL2)*gamma+L3*(gamma+delta-epsilon)-2*dL3-L3) / (q-alpha*z);
      
        err1 = min(abs(val-val2),abs(val-val3));

      else

        err1 = Inf;
      
      end

      if (abs(q-alpha*z)<0.01)||(err1<eps)
      
        err2 = abs(L1)*eps*N + abs(ckm0)*sqrt(numb-N+1) + abs(L2)*eps*(numb-N+1) + ...
               abs(log(z)) * ( abs(skm0)*sqrt(numb-N+1) + abs(L3)*eps*(numb-N+1) );
        err =  min(err1,err2);
      
      else
      
        err = err1;
      
      end
       
    end

  end

end


function co = findcoef4HeunCs(q,alpha,gamma,delta,epsilon)

  n = round(1-gamma);
  
  recur = @(k,ckm1,ckm2) ...
      (ckm1*(-q+(k-1)*(gamma-epsilon+delta+k-2)) + ...
       ckm2*((k-2)*epsilon+alpha))/(k*(gamma+k-1));
       
  ckm1 = 1; ckm2 = 0;
  
  co = epsilon^n/factorial(n);
    
  for k=1:n-1
    ckm0 = recur(k,ckm1,ckm2);
    co = co + ckm0 * epsilon^(n-k)/factorial(n-k);
    ckm2 = ckm1; ckm1 = ckm0; 
  end

end
