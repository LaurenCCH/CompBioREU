function[]=compare_plots_exp(num_sims,data_nums, data_cell,MLE_q_analytic)
%This function loops over ordered triples (data_nums,bw,scale_small_probs)and
%plots num_sims verusus abs(MLE_q_approx_simulation-MLE_q_analytic).
p_k=zeros(length(data_nums));
c=zeros(length(data_nums));
p_c=zeros(length(data_nums));


        for k=1:length(data_nums)
           
            [p_k(k),k_const,SSE_k]=convergence_rate(num_sims,data_cell{k}.avg_ML_error);
            
            [c(k),p_c(k),SSE_c] = convergence_rate_C(data_cell{k}.avg_ML_error,num_sims);

            
            set(figure,'DefaultFigureWindowStyle','docked')
            hold on
            plot(num_sims,data_cell{k}.avg_ML_error,'o')
            title(sprintf("Number of simulations versus average absolute error and variance (data_nums=%f)" , data_nums(k)), 'Interpreter', 'none');
            xlabel(('Num_sims') , 'Interpreter', 'none')
            ylabel('MLE error and MLE variance') 
            plot(num_sims,data_cell{k}.sample_error_variance,'o')
            fine_range_mesh=num_sims(end):(.01):num_sims(1);
            
            %This is the plot for convergence_rate
            plot(fine_range_mesh, (k_const./(fine_range_mesh).^p_k(k)),'r')
            
            plot(fine_range_mesh, (1./(fine_range_mesh).^p_c(k))+c(k),'b')
            %This annotation is for convergence_rate
            annotation('textbox','String',sprintf("p_k=%f\nk\\_const=%f\nSSE_k=%f\np=%f\nc=%f\nSSE_c=%f", p_k(k),k_const,SSE_k,p_c(k),c(k),SSE_c),'FitBoxToText','on');
            %annotation('textbox','String',sprintf("p=%f\nc=%f\nSSE=%f",p,c,SSE_c),'FitBoxToText','on');
            legend('avg\_ML\_error','sample\_error\_variance','decay fit k','decay fit c');
            
            hold off
            figure
            plot(num_sims,((data_cell{k}.avg_MLE_q_approx_simulation)-MLE_q_analytic(k)),'*')
            %ADD LABELS
            xlabel('Number of simulations');
            ylabel('Approximate MLE of q - MLE of q');
            annotation('textbox','String',sprintf("data\\_nums=%f", data_nums(k)),'FitBoxToText','on');
           
            hold off
           
            
            for num_sims_index=1:length(num_sims)    
                if num_sims_index==1 || num_sims_index==length(num_sims)||num_sims_index==ceil(length(num_sims)/2)
                        %sample_error_plot(data_cell{k}.e_sample_matrix(num_sims_index, :),data_nums(k), num_sims(num_sims_index));
                        %Generating sample error plot within this function
%                         histogram(data_cell{k}.e_sample_matrix(num_sims_index, :),20)
%                         annotation('textbox','String',sprintf("data_nums=%.2f\nnum_sims=%.2f", data_nums,num_sims),'FitBoxToText','on', 'Interpreter', 'none', 'Position', [.5, .7, .4, .2]);
%                         xlabel('Absolute error in the estimate of MLE q','FontSize', 16)
%                         ylabel('Frequency', 'FontSize', 16)
%                         title('Distribution of the absolute error in the estimate of the MLE of q over all samples', 'Interpreter', 'latex')
% 
%                         
                        hold off
               
                        figure
                        histogram((data_cell{k}.q_sample_matrix(num_sims_index, :)-MLE_q_analytic(k)), 10);
            
                        xlabel('Approximate MLE of q - MLE of q');
                        ylabel('Frequency')
                        title('Distribution of the error in the estimate of the MLE of q over all samples', 'Interpreter','latex')
                        annotation('textbox','String',sprintf("data\\_nums=%f", data_nums(k)), 'Position', [.5, .7, .4, .2],'FitBoxToText','on');
           
                end
            hold off
          
            
            end
            
            %histogram(data_cell{i,j,k}.e_sample_matrix,'BinWidth',0.015)
           
        end
    

hold off
figure

plot(data_nums, p_k, 'o')
xlabel('data_nums','Interpreter', 'none' );
ylabel('p_k')

figure
plot(data_nums, p_c, 'o')
xlabel('data_nums','Interpreter', 'none');
ylabel('p_c')

figure
plot(data_nums, c, 'o')
xlabel('data_nums', 'Interpreter', 'none');
ylabel('c')
end

