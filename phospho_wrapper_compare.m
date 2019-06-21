function [MLE_q_numeric,MLE_q_analytic,MLE_q_approx,numeric_LL,Max_LL,approx_LL,q_LL] = phospho_wrapper_compare(q,n)
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
Max_LL=likelihood(t,MLE_q_analytic);
q0 = q;
b = 0;
A = -1;
%We have the constraint A*q <= b, i.e. q>=0.
negLL=@(q)-1*likelihood(t,q);
options=optimset('MaxFunEvals',1000,'MaxIter',1000,'Display','iter');
[MLE_q_numeric,numeric_LL] = fmincon(negLL,q0,A,b,[],[],[],[],[],options);
numeric_LL=-numeric_LL;
num_sims=[100, 1000, 10000];
%num_sims=[100, 1000];
bw = [0.1,0.01];
scale_small_probs=[10,100];
MLE_q_approx=zeros(length(num_sims),2,2);
approx_LL=zeros(length(num_sims),2,2);
for i=1:length(num_sims)
    
    fprintf("%f percent complete", 100.0*i/length(num_sims));
    
    for j=1:2
        for k=1:2
            sth = 1/(num_sims(i)*scale_small_probs(k));           
            negLL_approx=@(q)-1*likelihood_approx(t,q,num_sims(i),bw(j),sth);
            [MLE_q_approx(i,j,k),approx_LL(i,j,k)] = fmincon(negLL_approx,q0,A,b,[],[],[],[],[],options);
            approx_LL(i,j,k)=-approx_LL(i,j,k);
        end
    end
end

q_LL=likelihood(t,q);

q_values=(MLE_q_analytic/2):0.01:(2*MLE_q_analytic);
plot(q_values,likelihood(t,q_values));
saveas(gcf,'likelihood_plot');
end