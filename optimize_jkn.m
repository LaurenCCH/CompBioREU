function [localstruct]=optimize_jkn(j,k,n_index,num_sims,num_samples,scale_small_probs,bw,q0,A,b,options,t,MLE_q_analytic)


% CREATE A STRUCT HERE
%data_cell{j,k,n_index} = struct();
localstruct = struct();


% ADD ALL OF THE APPROPRIATE ATTRIBUTE FIELDS
localstruct.avg_ML_error = zeros(size(num_sims));
localstruct.avg_MLE_q_approx_simulation = zeros(size(num_sims));
localstruct.avg_approx_LL_simulation = zeros(size(num_sims));
localstruct.sum_e_sample_squared = zeros(size(num_sims));


for i=1:length(num_sims)
    avg_ML_error = 0;
    avg_MLE_q_approx_simulation = 0;
    avg_approx_LL_simulation = 0;
    sum_e_sample_squared = 0;
 
    parfor num_samp_index=1:num_samples(i)

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
        e_sample=(abs(q_sample-(MLE_q_analytic(n_index))));

        % avg_ML_error(i,j,k,n_index)=avg_ML_error(i,j,k,n_index)+e_sample;
        avg_ML_error = avg_ML_error+e_sample;
        %Save e_samples here, make array, plot into a histogram outside loop.
        
        
        % avg_MLE_q_approx_simulation(i,j,k,n_index)=avg_MLE_q_approx_simulation(i,j,k,n_index)+q_sample;
        avg_MLE_q_approx_simulation = avg_MLE_q_approx_simulation+q_sample;

        % avg_approx_LL_simulation(i,j,k,n_index)=avg_approx_LL_simulation(i,j,k,n_index)+LL_sample;
        avg_approx_LL_simulation = avg_approx_LL_simulation+LL_sample;

        %sum_e_sample_squared(i,j,k,n_index)=sum_e_sample_squared(i,j,k,n_index) +((e_sample)^2);
        sum_e_sample_squared = sum_e_sample_squared+((e_sample)^2);

        % Remember that we have to flip the LL back becaus the
        % optimizer minimizes!

    end
    % avg_MLE_q_approx_simulation(i,j,k,n_index)=avg_MLE_q_approx_simulation(i,j,k,n_index)/num_samples(i);
    avg_MLE_q_approx_simulation = avg_MLE_q_approx_simulation/num_samples(i);

    % avg_approx_LL_simulation(i,j,k,n_index)=avg_approx_LL_simulation(i,j,k,n_index)/num_samples(i);
    avg_approx_LL_simulation =  avg_approx_LL_simulation/num_samples(i);

    % avg_ML_error(i,j,k,n_index)=avg_ML_error(i,j,k,n_index)/num_samples(i);
    avg_ML_error = avg_ML_error/num_samples(i);

    % not sure about this one
    % sample_error_variance(i,j,k,n_index)=(sum_e_sample_squared(i,j,k,n_index)/num_samples(i))-((avg_ML_error(i,j,k,n_index))^2);
    sample_error_variance = sum_e_sample_squared/num_samples(i)-((avg_ML_error)^2);


    % HERE WE ASSIGN THE STRUCT TO OUR CELL ARRAY
    localstruct.avg_ML_error(i) = avg_ML_error;
    localstruct.avg_MLE_q_approx_simulation(i) = avg_MLE_q_approx_simulation;
    localstruct.avg_approx_LL_simulation(i) = avg_approx_LL_simulation;
    localstruct.sum_e_sample_squared(i) = sum_e_sample_squared;
    localstruct.sample_error_variance(i) = sample_error_variance;
end
%histogram

end