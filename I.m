function[x_hats_to_sum]=I(x, x_i_hat)
prob_search=find (abs(x-x_i_hat)<.01*2*sqrt(3)*x_i_hat);
how_many_x_hats=size(prob_search);
x_hats_to_sum=((1/(.01*2*sqrt(3)*x_hat))*ones(how_many_x_hats));
end
