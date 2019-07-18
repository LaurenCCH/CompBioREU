

function[density] = Gaussian(sim_data, exp_data_i, cv)

bin_width = (cv/('\mu'))*sqrt(pi*2);
density = 0;

if abs(sim_data-exp_data_i)<bin_width
    density = (1/bin_width) * exp(1/2*((x-('\mu')/(cv/('\mu'))))^2); 
    
end 
end 
