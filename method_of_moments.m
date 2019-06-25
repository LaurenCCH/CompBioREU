%This code calcualtes the value of q for which the expected value of qe^-qt
%is equal to the sample mean.
function[q_mom]=method_of_moments(t)
sample_mean=mean(t);
%The expected value of qe^-qt is 1/q;
q_mom=1/sample_mean;
end