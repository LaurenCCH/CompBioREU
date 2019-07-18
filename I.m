function[density]=I(sim_data, exp_data_i, cv)
%evaluate the uniform kernel centered at x_hat at the data point x_i

%the bin width is chosen so that the ratio of the std dev to the mean is
%fixed at .01

bin_width=cv*2*sqrt(3)*sim_data;
density=0;
%if exp_data_i falls into the bin the density of the kernel at exp_data_i
%is the height of the bin, 1/bin_width. Otherwise, the density of the
%kernel at exp_data_i is 0.
if abs(sim_data-exp_data_i)<bin_width
    density=1/bin_width;
end
end
