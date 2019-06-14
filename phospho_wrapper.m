function [MLE_q,Max_LL,q_LL] = phospho_wrapper(q,n)
%length of small time step for simulating phosphorylation
timestep=.01;
[t]=phospho_times(q,timestep,n);

MLE_q_true=1/mean(t);
MLE_q_numeric=


Max_LL=likelihood(t,MLE_q_true);
q_LL=likelihood(t,q);

q_values=MLE_q/2:.01:2*MLE_q;
plot(q_values,likelihood(t,q_values));
saveas(gcf,'likelihood_plot');
end