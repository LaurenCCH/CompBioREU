function [MLE_q_numeric,MLE_q_analytic,Max_LL,q_LL] = phospho_wrapper(q,n)
%generate a set of synthetic data (t) for the first event (phosphorylation) 
%times for a Poisson process, given the provided rate, q, 
% and over a provided number of trials, n, and calculate the analytic MLE
% of q from that synthetic data and also calculate the numeric MLE
% of q and also return the log likelihood values of both estimators of q.
% MLE_q_numeric represents the numerical MLE of q
% MLE_q_analytic represents the analytic MLE of q
% Max_LL is the log likelihood of MLE_q_true
% q_LL is the log likelihood of q

%length of small time step for simulating phosphorylation
timestep=.01;
% this function returns the first phosphoralations times with n trials.
% q represets the rate of phosphoralations.
[t]=phospho_times(q,timestep,n);

MLE_q_analytic=1/mean(t);
q0 = q;
b = 0;
A = -1;
%We have the constraint A*q <= b, i.e. q>=0.
negLL=@(q)-1*likelihood(t,q);
options=optimset('MaxFunEvals',1000,'MaxIter',1000,'Display','iter');
[MLE_q_numeric] = fmincon(negLL,q0,A,b,[],[],[],[],[],options);

%,LL_numeric

Max_LL=likelihood(t,MLE_q_analytic);
q_LL=likelihood(t,q);

q_values=(MLE_q_analytic/2):0.01:(2*MLE_q_analytic);
plot(q_values,likelihood(t,q_values));
saveas(gcf,'likelihood_plot');
end