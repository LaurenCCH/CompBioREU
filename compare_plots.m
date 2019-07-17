function[]=compare_plots(num_sims,data_nums, bw, scale_small_probs,data_cell)
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
            title(sprintf("Num_sims versus average error and variance for\nbw=%f\nscale_small_probs=%f\ndata_nums=%f",bw(i) ,scale_small_probs(j) , data_nums(k)), 'Interpreter', 'none');
            xlabel(('Num_sims') , 'Interpreter', 'none')
            ylabel('MLE error and MLE variance') 
            plot(num_sims,data_cell{i,j,k}.sample_error_variance,'o')
            fine_range_mesh=num_sims(end):(.01):num_sims(1);
            
            %This is the plot for convergence_rate
            plot(fine_range_mesh, (k_const./(fine_range_mesh).^p_k(i,j,k)),'r')
            
            plot(fine_range_mesh, (1./(fine_range_mesh).^p_c(i,j,k))+c(i,j,k),'b')
            %This annotation is for convergence_rate
            annotation('textbox','String',sprintf("p_k=%f\nk\\_const=%f\nSSE_k=%f\np=%f\nc=%f\nSSE_c=%f", p_k(i,j,k),k_const,SSE_k,p_c(i,j,k),c(i,j,k),SSE_c),'FitBoxToText','on');
            %annotation('textbox','String',sprintf("p=%f\nc=%f\nSSE=%f",p,c,SSE_c),'FitBoxToText','on');
            legend('avg\_ML\_error','sample\_error\_variance','decay fit k','decay fit c');
            
            hold off
            figure
            
            for num_sims_index=1:length(num_sims)    
                if num_sims_index==1 || num_sims_index==length(num_sims)||num_sims_index==ceil(length(num_sims)/2)
                    sample_error_plot(data_cell{i,j,k}.e_sample_matrix(num_sims_index, :),data_nums(k),bw(i),scale_small_probs(j), num_sims(num_sims_index));
                end
            hold off
            figure
            
            end
            
            %histogram(data_cell{i,j,k}.e_sample_matrix,'BinWidth',0.015)
           
        end
    end
end
hold off
figure
compare_plots_3D(data_nums, scale_small_probs, bw, p_k, p_c, c);

          

end

