function[q_mom]=method_of_moments(t)
sample_mean=mean(t);
expectation=1/q;
q_mom=solve(expectation==sample_mean,q);
end