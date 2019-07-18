function[density] = Exponential(sim_data, exp_data_i, cv)

bin_width = cv*exp((-1/cv)*x);
density = 0;

if abs(sim_data-exp_data_i)<bin_width
    density = 1/bin_width;
    
end 
end 