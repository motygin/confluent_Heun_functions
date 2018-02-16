% confluent Heun function, a solution of the equation
% HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
%
% asymptotic expansion at z=infinity
% the first, power solution
%
% Usage:
% [val,dval,err,numb,wrnmsg] = HeunCinfA(q,alpha,gamma,delta,epsilon,z)
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
% Oleg V. Motygin, copyright 2017-2018, license: GNU GPL v3
%
% 20 December 2017
%
function [val,dval,err,numb,wrnmsg] = HeunCinfA(q,alpha,gamma,delta,epsilon,z)

  global Heun_asympt_klimit;
  
  if isempty(Heun_asympt_klimit)
    HeunOpts();
  end
  
  recurr1 = @(n,cnm1) cnm1*n/(z*epsilon)*(1+(-q+alpha/epsilon*(2*n-gamma-delta-1+alpha/epsilon)+ ...
      (gamma-epsilon+delta+1)*(1-n)+alpha-1)/n^2);

  recurr = @(n,cnm1,cnm2) recurr1(n,cnm1) + cnm2/(z^2*epsilon)* ...
      ((n-2+alpha/epsilon)*(gamma-n+1-alpha/epsilon))/n;

  wrnmsg = '';
  
  val = 1; dval = 0;
  err = 0; numb = 1;
  
  cnm3 = Inf; cnm2 = 1; cnm1 = recurr1(1,cnm2);
  dnm3 = Inf; dnm2 = 0; dnm1 = -cnm1/z;

  val = cnm2 + cnm1; dval = dnm1;
  
  vm3 = NaN; vm2 = NaN; vm1 = NaN; vm0 = val;
  dvm3 = NaN; dvm2 = NaN; dvm1 = NaN; dvm0 = dval;  

  numb = 2; small = sqrt(eps); 
  
  growcn = false; growdn = false;
  valstab = false; dvalstab = false;
  
  while (numb<=Heun_asympt_klimit) && ((abs(cnm3)>small) || ~(growcn||valstab) || ~(growdn||dvalstab))

    cnm0 = recurr(numb,cnm1,cnm2); dnm0 = -numb*cnm0/z;
    val = val + cnm0; dval = dval + dnm0;
    err = abs(cnm2);
    numb = numb+1;
    
    growcn = growcn || ((abs(cnm0)>abs(cnm1))&&(abs(cnm1)>abs(cnm2))&&(abs(cnm2)>abs(cnm3)));
    valstab = valstab || ((vm3==vm2)&&(vm2==vm1)&&(vm1==val));

    growdn = growdn || ((abs(dnm0)>abs(dnm1))&&(abs(dnm1)>abs(dnm2))&&(abs(dnm2)>abs(dnm3)));
    dvalstab = dvalstab || ((dvm3==dvm2)&&(dvm2==dvm1)&&(dvm1==dval));
    
    if ((abs(cnm2)>small) || ~(growcn||valstab))
      cnm3 = cnm2; cnm2 = cnm1; cnm1 = cnm0;
      vm3 = vm2; vm2 = vm1; vm1 = vm0; vm0 = val;
    end
  
    if ((abs(cnm2)>small) || ~(growdn||dvalstab))
      dnm3 = dnm2; dnm2 = dnm1; dnm1 = dnm0;
      dvm3 = dvm2; dvm2 = dvm1; dvm1 = dvm0; dvm0 = dval;
    end
    
  end  
  
  val = (-z)^(-alpha/epsilon) * vm3;
  dval = (-z)^(-alpha/epsilon) * (dvm3-alpha/epsilon*vm3/z);
  err = abs(z^(-alpha/epsilon)) * err;
  
  if ( isinf(val) || isinf(dval) || isnan(val) || isnan(dval) )
    
    wrnmsg = 'HeunCinfA: failed convergence of recurrence and summation; '; 
    val = NaN; dval = NaN; err = NaN;
      
  end

end

    
