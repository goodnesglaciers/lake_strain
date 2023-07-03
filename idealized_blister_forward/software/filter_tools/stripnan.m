function b = stripnan(a)
%b = stripnan(a)
%strips out all NaN values from vector a
%Only works if a is a vector
%If a is a column vector, stripnan returns a column vector
%If a is a row vector, stripnan returns a row vector

[n, m] = size(a);
if all([n, m] ~= 1)
   return
end

if (m == 1)
   j = 0;
   for i = 1:n
      if ~isnan(a(i,1))
         j = j + 1;
         b(j,1) = a(i,1);
      end
   end
else
   j = 0;
   for i = 1:m
      if ~isnan(a(1,i))
         j = j + 1;
         b(1,j) = a(1,i);
      end
   end
end
