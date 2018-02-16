% for confluent Heun function, a solution of the equation
% HeunC''(z)+(gamma/z+delta/(z-1)+epsilon)*HeunC'(z)+(alpha*z-q)/(z*(z-1))*HeunC(z) = 0
%
% HeunCjoin0infA finds connection coefficients C0, Cs, such that
% C0 * HeunC00(z) + Cs * HeunCs00(z) analytically continues to
% the first, power solution at infinity \exp(i\theta) \infty
% (see HeunCinfA)
%
% Usage:
% [C0,Cs,err,numb,wrnmsg] = HeunCjoin0infA(q,alpha,gamma,delta,epsilon,theta)
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
function [C0,Cs,err,numb,wrnmsg] = HeunCjoin0infA(q,alpha,gamma,delta,epsilon,theta,varargin)
  
  wrnmsg = '';
  
  [C0,Cs,err,numb,consts_known] = extrdatfromsav(q,alpha,gamma,delta,epsilon,theta);
  
  if consts_known
  
    numb = 0;
    
  else
  
    [R,N] = findR();
  
    R0 = R/(abs(eps)+abs(epsilon));
  
    infpt = 2 * R0 * exp(1j*theta);
    joinpt = 1j * min(1,R0) * sign(imag(exp(1j*theta)));
  
    [valinf,dvalinf,errinf,numbinf,wrnmsginf] = HeunCinfA(q,alpha,gamma,delta,epsilon,infpt);
    [valJinf,dvalJinf,errJinf,numbJinf,wrnmsgJinf] = HeunCconnect(q,alpha,gamma,delta,epsilon,joinpt,infpt,valinf,dvalinf);
    
    [valJ0,dvalJ0,errJ0,numbJ0,wrnmsgJ0] = HeunC0(q,alpha,gamma,delta,epsilon,joinpt,varargin);
    [valJs,dvalJs,errJs,numbJs,wrnmsgJs] = HeunCs0(q,alpha,gamma,delta,epsilon,joinpt);
  
    wrnmsg = strcat(wrnmsginf,wrnmsgJinf,wrnmsgJ0,wrnmsgJs);
    err = errinf + errJinf + errJ0 + errJs;
    numb = numbinf + numbJinf + numbJ0 + numbJs;
  
    if isempty(wrnmsg)
      m = [valJ0,valJs;dvalJ0,dvalJs];
      b = [valJinf;dvalJinf];
      cl = m \ b;
      C0 = cl(1); Cs = cl(2);    
    end

    keepdattosav(q,alpha,gamma,delta,epsilon,theta,C0,Cs,err,numb);

  end
  
end

function [C0,Cs,err,numb,consts_known] = extrdatfromsav(q,alpha,gamma,delta,epsilon,theta)

  global savdata0inf;

  C0 = NaN; Cs = NaN; err = NaN; numb = 0; consts_known = false;
  if length(savdata0inf)~=0
    for k=1:length(savdata0inf)
      if (savdata0inf(k).q==q)&&(savdata0inf(k).alpha==alpha)&&(savdata0inf(k).gamma==gamma)&& ...
         (savdata0inf(k).delta==delta)&&(savdata0inf(k).epsilon==epsilon)&&(savdata0inf(k).theta==theta)
        C0 = savdata0inf(k).C0;
        Cs = savdata0inf(k).Cs;
        err = savdata0inf(k).err;
        numb = savdata0inf(k).numb;
        consts_known = true;
        break;
      end
    end
  end
end

function keepdattosav(q,alpha,gamma,delta,epsilon,theta,C0,Cs,err,numb);
 
  global savdata0inf;
  
  global Heun_memlimit;
  
  if isempty(Heun_memlimit)
    HeunOpts();
  end

  if length(savdata0inf)==0
    savdata0inf=struct('q',q,'alpha',alpha,'gamma',gamma,'delta',delta,...
    'epsilon',epsilon,'theta',theta,'err',err,'numb',numb,'C0',C0,'Cs',Cs);
  else
    if length(savdata0inf)<=Heun_memlimit
      savdata0inf(end+1)=struct('q',q,'alpha',alpha,'gamma',gamma,'delta',delta,...
      'epsilon',epsilon,'theta',theta,'err',err,'numb',numb,'C0',C0,'Cs',Cs);
    else
      savdata0inf(1).q=q; savdata0inf(1).alpha=alpha;
      savdata0inf(1).gamma=gamma; savdata0inf(1).delta=delta;
      savdata0inf(1).epsilon=epsilon; savdata0inf(1).theta=theta;
      savdata0inf(1).err=err; savdata0inf(1).numb=numb;
      savdata0inf(1).C0=C0; savdata0inf(1).Cs=Cs;
      savdata0inf = shift(savdata0inf,-1);
    end  
  end
end
