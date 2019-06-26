function[]=compare_plots(avg_MLE_q_approx_simulation,num_sims,data_nums, bw, scale_small_probs,MLE_q_analytic, sample_error_variance)
%This function loops over ordered triples (data_nums,bw,scale_small_probs)and
%plots num_sims verusus abs(MLE_q_approx_simulation-MLE_q_analytic).

avg_ML_error=zeros(length(num_sims),length(bw),length(scale_small_probs),length(data_nums));
for i=1:length(bw)
    for j=1:length(scale_small_probs)
        for k=1:length(data_nums)
           
            avg_ML_error(:,i,j,k)=(abs(avg_MLE_q_approx_simulation(:,i,j,k)-(MLE_q_analytic)));
            set(figure,'DefaultFigureWindowStyle','docked')
            plot(num_sims,avg_ML_error(:,i,j,k),'o')
            string_num_sims=num2str(num_sims);
            title(("Num_sims versus error for num_sims="+string_num_sims+" , bw=" +bw(i)+" , scale_small_probs="+scale_small_probs(j)+" , and data_nums="+data_nums(k)), 'Interpreter', 'none')  
            xlabel(('Num_sims') , 'Interpreter', 'none')
            ylabel('ML_error') 
            
            hold on
            
            plot(num_sims,sample_error_variance(:,i,j,k),'o')
            string_num_sims=num2str(num_sims);
            title(("Num_sims versus error variance for num_sims="+string_num_sims+" , bw=" +bw(i)+" , scale_small_probs="+scale_small_probs(j)+" , and data_nums="+data_nums(k)), 'Interpreter', 'none')  
            xlabel(('Num_sims') , 'Interpreter', 'none')
            ylabel('Error variance') 
            
            
        end
    end
end
end

