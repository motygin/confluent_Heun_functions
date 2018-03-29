% Change internal parameters of HeunC, HeunCs, HeunC0, HeunCs0, HeunC00, HeunCs00, 
% HeunCfromZ0, HeunCconnect, HeunCinfA, HeunCjoin0infA, HeunCjoin1infA, HeunCjoin01
%
% Usage: HeunOpts('cont_coef',Heun_cont_coef,'klimit',Heun_klimit,'optserterms',Heun_optserterms,
%  'asympt_klimit',Heun_asympt_klimit,'proxco',Heun_proxco,'proxcoinfrel',Heun_proxcoinf_rel,'memlimit',Heun_memlimit)
%
% parameters are optional
%
% if a pair of parameters is omitted, but the corresponding global variable
% is empty then it is set to the default value
%
% Use [] as the value to reset coefficient to its default
% e.g. HeunOpts('cont_coef',[])
%
% Heun_cont_coef: for each power expansion of the analytic continuation procedure 
% the coefficient is the relative (to the radius of convergence) distance from 
% the centre to the calculated point.
% By default Heun_cont_coef = 0.38
%
% Heun_klimit is the maximum number of power series' terms.
% Default value is 1000 
%
% Heun_optserterms is the number of power series' terms considered as in some sense optimal.
% Default value is 40
%
% Heun_asympt_klimit is the maximum number of asymptotic series' terms.
% Default value is 200
%
% Heun_proxco and Heun_proxcoinf_rel specify relative proximity
% to singular point where special representation is used
% By default Heun_proxco = 0.05; Heun_proxcoinf_rel = 1.0 
%
% Heun_memlimit is the maximum number of sets of data (parameters of confluent Heun function 
% and corresponding connection coefficients) which are kept in memory.
% Default value is 500
%
% Oleg V. Motygin, copyright 2018, license: GNU GPL v3
%
% 27 March 2018
%
function HeunOpts(varargin)

  global Heun_cont_coef Heun_klimit Heun_optserterms Heun_asympt_klimit Heun_proxco Heun_proxcoinf_rel Heun_memlimit;
  
  [reg, props] = parseparams(varargin);
  opts = cell2struct(props(2:2:end),props(1:2:end),2);
  
  if isfield(opts,'cont_coef')
    if isempty(opts.cont_coef)
      Heun_cont_coef = 0.38;
    else
      Heun_cont_coef = opts.cont_coef;
    end
  else
    if isempty(Heun_cont_coef)
      Heun_cont_coef = 0.38;
    end
  end

  if isfield(opts,'klimit')
    if isempty(opts.klimit)
      Heun_klimit = 1000;
    else
      Heun_klimit = opts.klimit;
    end
  else
    if isempty(Heun_klimit)
      Heun_klimit = 1000;
    end
  end
  
  if isfield(opts,'optserterms')
    if isempty(opts.optserterms)
      Heun_optserterms = 40;
    else
      Heun_optserterms = opts.optserterms;
    end
  else
    if isempty(Heun_optserterms)
      Heun_optserterms = 40;
    end
  end
  
  if isfield(opts,'asympt_klimit')
    if isempty(opts.asympt_klimit)
      Heun_asympt_klimit = 200;
    else
      Heun_asympt_klimit = opts.asympt_klimit;
    end
  else
    if isempty(Heun_asympt_klimit)
      Heun_asympt_klimit = 200;
    end
  end
  
  if isfield(opts,'proxco')
    if isempty(opts.proxco)
      Heun_proxco = 0.05;
    else
      Heun_proxco = opts.proxco;
    end
  else
    if isempty(Heun_proxco)
      Heun_proxco = 0.05;
    end
  end
  
  if isfield(opts,'proxcoinfrel')
    if isempty(opts.proxcoinfrel)
      Heun_proxcoinf_rel = 1.0;
    else
      Heun_proxcoinf_rel = opts.proxcoinfrel;
    end
  else
    if isempty(Heun_proxcoinf_rel)
      Heun_proxcoinf_rel = 1.0;
    end
  end
  
  if isfield(opts,'memlimit')
    if isempty(opts.memlimit)
      Heun_memlimit = 500;
    else
      Heun_memlimit = opts.memlimit;
    end
  else
    if isempty(Heun_memlimit)
      Heun_memlimit = 500;
    end
  end

end
