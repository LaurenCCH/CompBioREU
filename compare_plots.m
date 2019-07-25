function[]=compare_plots(num_sims,data_nums, bw, scale_small_probs,data_cell,MLE_q_analytic,num_samples)
%This function loops over ordered triples (data_nums,bw,scale_small_probs)and
%plots num_sims verusus abs(MLE_q_approx_simulation-MLE_q_analytic).
p_k=zeros(length(bw),length(scale_small_probs),length(data_nums));
c=zeros(length(bw),length(scale_small_probs),length(data_nums));
p_c=zeros(length(bw),length(scale_small_probs),length(data_nums));

for i=1:length(bw)
   for j=1:length(scale_small_probs)
        for k=1:length(data_nums)
           
            [p_k(i,j,k),k_const,SSE_k]=convergence_rate(num_sims,data_cell{i,j,k}.avg_ML_error);
            
            [c(i,j,k),p_c(i,j,k),SSE_c] = convergence_rate_C(data_cell{i,j,k}.avg_ML_error,num_sims);

              
            
            

            set(figure,'DefaultFigureWindowStyle','docked')
            hold on
            plot(num_sims,data_cell{i,j,k}.avg_ML_error,'o')
            title(sprintf("Number of simulations versus the average asolute error and its variance\n data size=%.2f\nbin width=%.2f\ndegenerate probability=(number of simulations)^{-1}%.2f\nnumber of samples =%.0f", data_nums(k),bw(i), 1/scale_small_probs(j),num_samples));
            xlabel(('Number of simulations') , 'Interpreter', 'none') 
            plot(num_sims,data_cell{i,j,k}.sample_error_variance,'o')
            fine_range_mesh=num_sims(end):(.01):num_sims(1);

            %This is the plot for convergence_rate
            plot(fine_range_mesh, (k_const./(fine_range_mesh).^p_k(i,j,k)),'r')

            plot(fine_range_mesh, (1./(fine_range_mesh).^p_c(i,j,k))+c(i,j,k),'b')
            %This annotation is for convergence_rate
            annotation('textbox','String',sprintf("p_k=%.2f\nk\\_const=%.2f\nSSE_k=%.2f\np=%.2f\nc=%.2f\nSSE_c=%.2f", p_k(i,j,k),k_const,SSE_k,p_c(i,j,k),c(i,j,k),SSE_c),'Position', [.3, .5, .4, .2],'FitBoxToText','on');
            %annotation('textbox','String',sprintf("p=%f\nc=%f\nSSE=%f",p,c,SSE_c),'FitBoxToText','on');
            legend('Average absolute error','Variance of the absolute error','decay fit k','decay fit c');
     

            hold off
            figure
            plot(num_sims,((data_cell{i,j,k}.avg_MLE_q_approx_simulation)-MLE_q_analytic(k)),'*')
            %ADD LABELS
            xlabel('Number of simulations');
            ylabel('average AMLE - MLE');
            annotation('textbox','String',sprintf("data size=%.0f\nbin width=%.2f\ndegenerate probability=(number of simulations)^{-1}%.2f\nnumber of samples =%.0f", data_nums(k),bw(i), 1/scale_small_probs(j),num_samples),'Position', [.3, .7, .4, .2],'FitBoxToText','on');
           
            hold off
           
            
            for num_sims_index=1:length(num_sims)    
                if num_sims_index==1 || num_sims_index==length(num_sims)||num_sims_index==ceil(length(num_sims)/2)
                        sample_error_plot(data_cell{i,j,k}.e_sample_matrix(num_sims_index, :),data_nums(k),bw(i),scale_small_probs(j), num_sims(num_sims_index));
                        
                        
                        
                        hold off
               
                        figure
                        histogram((data_cell{i,j,k}.q_sample_matrix(num_sims_index, :)-MLE_q_analytic(k)), 10);
            
                        xlabel('Approximate MLE of \lambda - MLE of \lambda');
                        ylabel('Frequency')
                        annotation('textbox','String',sprintf("data size=%.0f\nbin width=%.0f\ndegenerate probability=(number of simulations)^{-1}%.2f\nnumber of samples =%.0f", data_nums(k),bw(i), 1/scale_small_probs(j),num_samples), 'Position', [.5, .7, .4, .2],'FitBoxToText','on');
           
                end
            hold off
          
            
            end
            
%             %histogram(data_cell{i,j,k}.e_sample_matrix,'BinWidth',0.015)
           
        end
    end
end
hold off
figure
compare_plots_3D(data_nums, scale_small_probs, bw, p_k, p_c, c);

          

end

