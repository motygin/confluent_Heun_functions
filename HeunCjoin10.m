% for confluent Heun function, a solution of the equation
% HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
%
% HeunCjoin10 finds matrix C of connection coefficients, such that
% C(1,1) * HeunC1(z) + C(1,2) * HeunCs1(z) analytically continues to
% the first local solution at z=0 (HeunC0), and
% C(2,1) * HeunC1(z) + C(2,2) * HeunCs1(z) continues to the second
% local solution at z=0 (HeunCs0)
%
% Usage:
% [C10,err,numb,wrnmsg] = HeunCjoin10(q,alpha,gamma,delta,epsilon)
%
% Returned parameters:
% C10 is the matrix of connection coefficients
% err is the estimated error
% numb is the number of power series terms needed for the evaluation
% wrnmsg is a warning message:
%   it is empty if computations are ok
%   otherwise it is a diagnostic message and the function 
%   returns C10 = [NaN NaN; NaN NaN]
%
% Oleg V. Motygin, copyright 2018, license: GNU GPL v3
%
% 15 March 2018
%
function [C10,err,numb,wrnmsg] = HeunCjoin10(q,alpha,gamma,delta,epsilon)
    
  wrnmsg = '';
  
  [C10,err,numb,consts_known] = extrdatfromsav(q,alpha,gamma,delta,epsilon);
  
  if consts_known
  
    numb = 0;
    
  else

    joinpt = 0.5;
  
    [val0,dval0,err0,numb0,wrnmsg0] = HeunC0(q,alpha,gamma,delta,epsilon,joinpt);
    [val0s,dval0s,err0s,numb0s,wrnmsg0s] = HeunCs0(q,alpha,gamma,delta,epsilon,joinpt);
    [val1,dval1,err1,numb1,wrnmsg1] = HeunC1(q,alpha,gamma,delta,epsilon,joinpt);
    [val1s,dval1s,err1s,numb1s,wrnmsg1s] = HeunCs1(q,alpha,gamma,delta,epsilon,joinpt);

    wrnmsg = strcat(wrnmsg0,wrnmsg0s,wrnmsg1,wrnmsg1s);
    err = err0/(1+abs(val0)) + err0s/(1+abs(val0s)) + err1/(1+abs(val1)) + err1s/(1+abs(val1s));
    numb = numb0 + numb0s + numb1 + numb1s;
  
    if isempty(wrnmsg)
      m = [val1,val1s;dval1,dval1s];
      b = [val0 val0s;dval0 dval0s];
      C10 = transpose(m \ b);
      
      err = err + sum(sum(abs(m*transpose(C10)-b)));
      
    end

    keepdattosav(q,alpha,gamma,delta,epsilon,C10,err,numb);

  end
  
end

function [val,dval,err,numb,wrnmsg] = HeunC1(q,alpha,gamma,delta,epsilon,z)
    [val,dval,err,numb,wrnmsg] = HeunC0(q-alpha,-alpha,delta,gamma,-epsilon,1-z);
    dval = -dval;
end
   
function [val,dval,err,numb,wrnmsg] = HeunCs1(q,alpha,gamma,delta,epsilon,z)
    [val,dval,err,numb,wrnmsg] = HeunCs0(q-alpha,-alpha,delta,gamma,-epsilon,1-z);
    dval = -dval;
end

function [C10,err,numb,consts_known] = extrdatfromsav(q,alpha,gamma,delta,epsilon)

  global savdata10;

  C10 = [NaN NaN; NaN NaN]; err = NaN; numb = 0; consts_known = false;
  if length(savdata10)~=0
    for k=1:length(savdata10)
      if (savdata10(k).q==q)&&(savdata10(k).alpha==alpha)&&(savdata10(k).gamma==gamma)&& ...
         (savdata10(k).delta==delta)&&(savdata10(k).epsilon==epsilon)
        C10 = savdata10(k).C10;
        err = savdata10(k).err;
        numb = savdata10(k).numb;
        consts_known = true;
        break;
      end
    end
  end
end

function keepdattosav(q,alpha,gamma,delta,epsilon,C10,err,numb);
 
  global savdata10;

  global Heun_memlimit;
  
  if isempty(Heun_memlimit)
    HeunOpts();
  end

  if length(savdata10)==0
    savdata10=struct('q',q,'alpha',alpha,'gamma',gamma,'delta',delta,...
    'epsilon',epsilon,'err',err,'numb',numb,'C10',C10);
  else
    if length(savdata10)<=Heun_memlimit
      savdata10(end+1)=struct('q',q,'alpha',alpha,'gamma',gamma,'delta',delta,...
      'epsilon',epsilon,'err',err,'numb',numb,'C10',C10);
    else
      savdata10(1).q=q; savdata10(1).alpha=alpha;
      savdata10(1).gamma=gamma; savdata10(1).delta=delta;
      savdata10(1).epsilon=epsilon;
      savdata10(1).err=err; savdata10(1).numb=numb;
      savdata10(1).C10=C10;
      savdata10 = shift(savdata10,-1);
    end  
  end
end
