% in HeunC*, HeunCs*, it is assumed that for |z|>2*R, asymptotic series 
% of type \sum_{n=0}^{\infty} b_n n!/z^n can be computed as superasymptotic
% with accuracy better than machine epsilon

function [R,N] = findR()

  global saveR saveN;

  if ~isempty(saveR) && ~isempty(saveN)

    R = saveR;
    N = saveN;
  
  else

    logeps = log(eps);
    R = -logeps;
    fact = 1;
    n = 1;
  
    while (true)
      n = n + 1;
      fact = fact * n;
      R0 = R;
      R = (log(fact)-logeps)/n;
      if (R > R0)
        break;
      end
    end
  
    N = n-1;
    R = (fact/n/eps)^(1/N);
    saveR = R; saveN = N;

  end
  
end
