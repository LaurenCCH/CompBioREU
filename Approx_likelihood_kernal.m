function[x_prob] = Approx_likelihood_kernal(x, num_sims, q_hat)
h = 1/(q_hat*100);
[x_hat] = phospho_times(q_hat,h,num_sims);
x_hat_sum=zeros(size(x));
for hat_index=1:length(x_hat)
    [x_hats_to_sum] = I(x,x_hat(hat_index));
    x_hat_scalar = log(sum(x_hats_to_sum)/num_sims);
    x_hat_sum(hat_index) = x_hat_scalar;
end
x_prob = sum(x_hat_sum);

% [I_function] = I(x,x_hat);
% t_prob(j) = 1/num_sims(sum(x_hat));
% sum(log(t_prob(j))) = sum(log(t_prob(j)));
end