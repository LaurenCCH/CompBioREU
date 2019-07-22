function[LL] = approx_likelihood_kernel(data, num_sims, q_hat,degenerate_prob, bw,kernel_type)
%This code approximates the log likelihood of given data (data).  
%Specifically, we use simulated data (sim_data) to generate an approximate 
%pdf based on the kernel provided by I. We calculate the LL of the given
%data provided the approximate pdf.
if strcmp(kernel_type,'exp')==1
    %the final argument represents cv and is unused for this kernel
    kernel_function=@(sim_data_j,exp_data_i,cv)Exponential(sim_data_j, exp_data_i,cv);
end
if strcmp(kernel_type,'indicator')==1
    %the final argument represents cv and is unused for this kernel
    kernel_function=@(sim_data_j,exp_data_i,cv)I(sim_data_j, exp_data_i,cv);
end
if strcmp(kernel_type,'gaussian')==1
    %the final argument represents cv and is unused for this kernel
    kernel_function=@(sim_data_j,exp_data_i,cv)Gaussian(sim_data_j, exp_data_i,cv);
end
if strcmp(kernel_type,'fixed_width_centered_window')==1
    %the final argument represents cv and is unused for this kernel
    kernel_function=@(sim_data_j,exp_data_i,cv)fixed_width_centered_window_kernel(sim_data_j, exp_data_i,bw);
end


%For simulating the continous stochastic process we choose the time step to
%be small relative to the occurrence rate, q.
h = 1/(q_hat*100);
cv=1/bw;
%simulated data for generating the approximate pdf
[sim_data] = phospho_times(q_hat,h,num_sims);
%the approximate pdf value of each point in data is stored in pdf_value
pdf_value=zeros(1,length(data));
%loop through data values to get corresponding approximate pdf values


%figure;
%hold on;


for data_index=1:length(data)
    %loop through the kernels centered on the simulated data to find the
    %density each kernel contributes to pdf_value(data_index). These values
    %sum to give the final pdf_value(data_index).
    for sim_data_index=1:num_sims
        pdf_value(data_index) = pdf_value(data_index)+kernel_function(sim_data(sim_data_index), data(data_index), cv);
    end
    pdf_value(data_index)=pdf_value(data_index)/(num_sims);
    if pdf_value(data_index)==0
        pdf_value(data_index)=degenerate_prob;
    end
    
    
    %plot(data(data_index),pdf_value(data_index),'*');
    
    
    
end
% title('Likelihood scatter')
% x_grid=0:0.01:max(data);
% plot(x_grid, Exponential(2, x_grid,42));
% 
% num_sims
% 
% for sim_data_index=1:num_sims
%     
%     this_kernel = zeros(1,length(x_grid));
%     for x_grid_point_idx=1:length(x_grid)
%        x_grid_point = x_grid(x_grid_point_idx);
%        this_kernel(x_grid_point_idx) = kernel_function(sim_data(sim_data_index),x_grid_point,cv);  
%     end
%     
%     plot(x_grid, this_kernel);
% end

LL=sum(log(pdf_value));
end