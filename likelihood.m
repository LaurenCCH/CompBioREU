% This function returns the loglikelihood of the parameter q given the data t.
% Here t is a vector of times for the first event of a process.
% We assume the process is a poisson process. 
% Under these assumptions t is known to follow an exponential distribution

% q is a row vector of rates
% t is a row vector of first event times
function [l] = likelihood(t,q)

% p(t,q) represents the probability distribution of t
p = @(t,q)(diag(q) * exp(-q' * t));

% l is what we are trying to return in this function which is the
% loglikelihood by summing the log values of p. 
l=sum(log(p(t,q)),2);

end





