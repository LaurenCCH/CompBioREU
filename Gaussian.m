

function[density] = Gaussian(sim_data_j, exp_data_i, cv)

sigma = (cv/sim_data_j);
density = (1/sigma*sqrt(pi*2)) * exp(-((exp_data_i-sim_data_j)^2)/(2*sigma^2)); 

end 
