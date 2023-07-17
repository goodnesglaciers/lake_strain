function [npwr2,n] = pwr2(x);
%
% If x is a vector PWR2 calculates the smallest integer 
% npwr2 s.t. 2^npwr2 >= length(x), and returns the value 
% npwr2 and 2^npwr2.
% If x is a scalar PWR2 calculates the smallest integer 
% npwr2 s.t. 2^npwr2 >= x, and returns the value npwr2
% and 2^npwr2.
% USAGE:  [npwr2,n] = pwr2(x)
%
%-------------------------------- jac

if (length(x) > 1)
    npts = length(x);
else
    npts = x;
end
n = 0;
while (1)
  if (2^n >= npts)
    break
  else
    n = n + 1;
  end
end
npwr2 = n;
n = 2^npwr2;
