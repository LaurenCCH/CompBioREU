function[]=compare_plots(avg_ML_error,sample_error_variance,num_sims,data_nums, bw, scale_small_probs)
%This function loops over ordered triples (data_nums,bw,scale_small_probs)and
%plots num_sims verusus abs(MLE_q_approx_simulation-MLE_q_analytic).

for i=1:length(bw)
    for j=1:length(scale_small_probs)
        for k=1:length(data_nums)
            [p,k_const,SSE]=convergence_rate(num_sims,avg_ML_error(:,i,j,k));
            set(figure,'DefaultFigureWindowStyle','docked')
            hold on
            plot(num_sims,avg_ML_error(:,i,j,k),'o')
            title(sprintf("Num_sims versus average error and variance for\nbw=%f\nscale_small_probs=%f\ndata_nums=%f",bw(i) ,scale_small_probs(j) , data_nums(k)), 'Interpreter', 'none');
            xlabel(('Num_sims') , 'Interpreter', 'none')
            ylabel('MLE error and MLE variance') 
            plot(num_sims,sample_error_variance(:,i,j,k),'o')
            fine_range_mesh=num_sims(1):0.01:num_sims(end);
            plot(fine_range_mesh, (k_const./fine_range_mesh).^p)
            annotation('textbox','String',sprintf("p=%f\nk\\_const=%f\nSSE=%f",p,k_const,SSE),'FitBoxToText','on');
            legend('avg\_ML\_error','sample\_error\_variance','decay fit');
        end
    end
end
end

