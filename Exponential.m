function[density] = Exponential(sim_data_j, exp_data_i,cv)

Gamma = (1/sim_data_j);
density = Gamma * exp(-Gamma * exp_data_i);


end 