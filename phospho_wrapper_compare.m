function [MLE_q_numeric,MLE_q_analytic,mom_q,numeric_LL,Max_LL,q_LL,mom_LL,data_nums,bw,scale_small_probs,num_sims,data_cell] = phospho_wrapper_compare(q)
%generate a set of synthetic data (t) for the first event (phosphorylation) 
%times for a Poisson process, given the provided rate, q, 
% and over a provided number of trials, n, and calculate the analytic MLE
% of q from that synthetic data and also calculate the numeric MLE
% of q and also return the log likelihood values of both estimators of q.
% MLE_q_numeric represents the numerical MLE of q
% MLE_q_analytic represents the analytic MLE of q
% Max_LL is the log likelihood of MLE_q_true
% q_LL is the log likelihood of q

gofast_mode=0;

% Make sure that the timestep, called h in other places, is always
% infintestimal compared to our provided value for q
timestep=1/(q*100);

if(gofast_mode==1)
    
    % These represent the number of samples, the number of t values we
    % generate
    data_nums = [100,1000];
    
    % These represetn the number of simulations that we will run to
    % approximate the pdf.
    num_sims=(100:100:300);
    
    % These represent the bin widths for approximating the pdf as a
    % histogram. The binwidth should be larger than the
    % inintesimal timestep used to simulate the process.
    bw = [timestep*10, timestep*100];
    
    % These scale our q values
    scale_small_probs=[10, 100];
    
    %num_samples=ceil(10000./num_sims);
    num_samples=10*ones(size(num_sims));
    
else
    % This is the go slow, full mode
    data_nums = [1000,500, 250,100, 50];
    num_sims=[50 250 500 750 1000 2500 5000 7500 10000];
    num_sims=fliplr(num_sims);
    bw = [timestep*2, timestep*5, timestep*10,timestep*20, timestep*30, timestep*40,timestep*50, timestep*60, timestep*70, timestep*80, timestep*90, timestep*100];
    scale_small_probs=[10,20, 30, 40, 50, 60, 70, 80, 90 100];
    num_samples=10*ones(size(num_sims));
end


data_cell = cell(length(bw),length(scale_small_probs),length(data_nums));

MLE_q_numeric=zeros(size(data_nums));
MLE_q_analytic=zeros(size(data_nums));
mom_q=zeros(size(data_nums));
numeric_LL=zeros(size(data_nums));
Max_LL=zeros(size(data_nums));
q_LL=zeros(size(data_nums));
mom_LL=zeros(size(data_nums));

for n_index=1:length(data_nums)
    n=data_nums(n_index);
 
    % this function returns the first phosphoralations times, t,  with n trials.
    % q represets the rate of phosphoralations.
    [t]=phospho_times(q,timestep,n);
   
    % Find the analytic MLE and best loglikelihoods for the provided t.
    MLE_q_analytic(n_index)=1/mean(t);
    Max_LL(n_index)=likelihood(t,MLE_q_analytic(n_index));
   
    % We want to find the maximum loglikelihood, but our optimizer wants to
    % minimize so we wrap the function to provide our objective function.
    negLL=@(q)-1*likelihood(t,q);
    
    % These are the options for the optimizer
    % We have the constraint A*q <= b, i.e. q>=0.
    q0 = q;
    b = 0;
    A = -1;
    options=optimset('MaxFunEvals',1000,'MaxIter',1000,'Display','iter');
    
    % Run the optimizer and find the numeric MLE and its associated LL.
    [MLE_q_numeric(n_index),numeric_LL(n_index)] = fmincon(negLL,q0,A,b,[],[],[],[],[],options);
    
    % Recall that our objective function was flipped, so flip the result
    % back.
    numeric_LL(n_index)=-numeric_LL(n_index);
    
    % The loglikelihood of the true parameter q
    q_LL(n_index)=likelihood(t,q);

    %the estimate of q by the method of moments
    mom_q(n_index)=method_of_moments(t);
    %get the loglikelihood of mom_q
    mom_LL(n_index)=likelihood(t,mom_q(n_index));

    % Prepare the domain values, q_values for which we will plot the
    % loglikelihoods
    q_values=(MLE_q_analytic(n_index)/2):0.01:(2*MLE_q_analytic(n_index));

    % Plot and save
    %compare_plots(MLE_q_approx_simulation,num_sims,data_nums, bw, scale_small_probs,MLE_q_analytic);
    plot(q_values,likelihood(t,q_values));
    %saveas(gcf,'likelihood_plot');
    
    % Loop through our desired number of simulations, indexed by i,
    % num_sims(i)
    
    % Loop through our desired bin widths, indexed by j, bw(j)
    
    
    
    
    for j=1:length(bw)
        % Loop through our desired small probability scaling factors,
        % indexed by k, scale_small_probs(k)
        for k=1:length(scale_small_probs)
            [localstruct]=optimize_jkn(j,k,n_index,num_sims,num_samples,scale_small_probs,bw,q0,A,b,options,t,MLE_q_analytic);
            data_cell{j,k,n_index} = localstruct;     
        end
    end
end

end
 