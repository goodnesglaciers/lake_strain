function factorial  = fact(n)
%

if n == 0
      factorial = 1;
   else
       i = 1:n;
       factorial = prod(i);
   end

