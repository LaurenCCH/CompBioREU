function [MLE_q_numeric,MLE_q_analytic,avg_MLE_q_approx_simulation,mom_q,numeric_LL,Max_LL,avg_approx_LL_simulation,q_LL,mom_LL,avg_ML_error,data_nums,bw,scale_small_probs,num_sims, sample_error_variance] = phospho_wrapper_compare(q)
%generate a set of synthetic data (t) for the first event (phosphorylation) 
%times for a Poisson process, given the provided rate, q, 
% and over a provided number of trials, n, and calculate the analytic MLE
% of q from that synthetic data and also calculate the numeric MLE
% of q and also return the log likelihood values of both estimators of q.
% MLE_q_numeric represents the numerical MLE of q
% MLE_q_analytic represents the analytic MLE of q
% Max_LL is the log likelihood of MLE_q_true
% q_LL is the log likelihood of q

gofast_mode=1;

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
    data_nums = [100,250,500, 750, 1000, 2500, 5000, 7500, 10000];
    num_sims=(100:100:10000);
    bw = [timestep*10, timestep*100];
    scale_small_probs=[10, 25, 50, 75, 100];
    num_samples=100*ones(size(num_sims));
end

% Preallocate the matrices holding the MLE and associated LL for each
% simulation run.
avg_MLE_q_approx_simulation=zeros(length(num_sims),length(bw),length(scale_small_probs),length(data_nums));
avg_approx_LL_simulation=zeros(length(num_sims),length(bw),length(scale_small_probs),length(data_nums));
avg_ML_error=zeros(length(num_sims),length(bw),length(scale_small_probs),length(data_nums));
sum_e_sample_squared=zeros(length(num_sims),length(bw),length(scale_small_probs),length(data_nums));
sample_error_variance=zeros(length(num_sims),length(bw),length(scale_small_probs),length(data_nums));
% Loop through our desired number of samples, indexed by n_index,
% data_nums(n_index)
for n_index=1:length(data_nums)
    n=data_nums(n_index);
 
    % this function returns the first phosphoralations times, t,  with n trials.
    % q represets the rate of phosphoralations.
    [t]=phospho_times(q,timestep,n);
   
    % Find the analytic MLE and best loglikelihoods for the provided t.
    MLE_q_analytic=1/mean(t);
    Max_LL=likelihood(t,MLE_q_analytic);
   
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
    [MLE_q_numeric,numeric_LL] = fmincon(negLL,q0,A,b,[],[],[],[],[],options);
    
    % Recall that our objective function was flipped, so flip the result
    % back.
    numeric_LL=-numeric_LL;
    
    
    
    % Loop through our desired number of simulations, indexed by i,
    % num_sims(i)
    
   

    % Loop through our desired bin widths, indexed by j, bw(j)
    for j=1:length(bw)

        % Loop through our desired small probability scaling factors,
        % indexed by k, scale_small_probs(k)
        for k=1:length(scale_small_probs)

             for i=1:length(num_sims)
                for num_samp_index=1:num_samples(i)

                    % Because we can underflow the probability, we need to
                    % sometimes replace the zero value with some degenerate
                    % nonzero probability to prevent the loglikeilhood from
                    % exploding. This value is scaled according to the number
                    % of simulations and the small probabillity scaling factor.
                    degenerate_probability = 1/(num_sims(i)*scale_small_probs(k));

                    % Prepare our objective function for optimization
                    negLL_approx=@(q)-1*likelihood_approx(t,q,num_sims(i),bw(j),degenerate_probability);

                    % Do the optimization and return the simulation MLE and LL
                    [q_sample,LL_sample] = fmincon(negLL_approx,q0,A,b,[],[],[],[],[],options);
                    LL_sample=-LL_sample;
                    e_sample=(abs(q_sample-(MLE_q_analytic)));
                    avg_ML_error(i,j,k,n_index)=avg_ML_error(i,j,k,n_index)+e_sample;

                    avg_MLE_q_approx_simulation(i,j,k,n_index)=avg_MLE_q_approx_simulation(i,j,k,n_index)+q_sample;
                    avg_approx_LL_simulation(i,j,k,n_index)=avg_approx_LL_simulation(i,j,k,n_index)+LL_sample;

                    sum_e_sample_squared(i,j,k,n_index)=sum_e_sample_squared(i,j,k,n_index) +((e_sample)^2);
                    % Remember that we have to flip the LL back becaus the
                    % optimizer minimizes!

                end
                avg_MLE_q_approx_simulation(i,j,k,n_index)=avg_MLE_q_approx_simulation(i,j,k,n_index)/num_samples(i);
                avg_approx_LL_simulation(i,j,k,n_index)=avg_approx_LL_simulation(i,j,k,n_index)/num_samples(i);
                avg_ML_error(i,j,k,n_index)=avg_ML_error(i,j,k,n_index)/num_samples(i);

                sample_error_variance(i,j,k,n_index)=(sum_e_sample_squared(i,j,k,n_index)/num_samples(i))-((avg_ML_error(i,j,k,n_index))^2);
             end
        end
    end
end

% The loglikelihood of the true parameter q
q_LL=likelihood(t,q);

%the estimate of q by the method of moments
mom_q=method_of_moments(t);
%get the loglikelihood of mom_q
mom_LL=likelihood(t,mom_q);

% Prepare the domain values, q_values for which we will plot the
% loglikelihoods
q_values=(MLE_q_analytic/2):0.01:(2*MLE_q_analytic);

% Plot and save
%compare_plots(MLE_q_approx_simulation,num_sims,data_nums, bw, scale_small_probs,MLE_q_analytic);
plot(q_values,likelihood(t,q_values));
%saveas(gcf,'likelihood_plot');
end
 