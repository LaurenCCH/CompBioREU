
function [l] = likelihood(t,q)
%q is a row vector of rates
%t is a row vector of times
p = @(t,q)(diag(q) * exp(-q' * t));

% this for loop calculates the likelehood 
% of the data
%L = 1
%for i=1 : n
%    L = L * p(t,q)
%    return L
%end


%t=0
%for i=1 : n
%   t = t + i
%   return t
%end

%q = @(t)(1/t)
fprintf("q is %f\n",q);
l=sum(log(p(t,q)),2);
size(l);








end
    
    
