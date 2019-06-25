function[]=compare_plots(MLE_q_approx_simulation,num_sims,data_nums, bw, scale_small_probs,MLE_q_analytic)
%This function loops over ordered triples (data_nums,bw,scale_small_probs)and
%plots num_sims verusus abs(MLE_q_approx_simulation-MLE_q_analytic).
for i=1:length(num_sims)
    for j=1:length(bw)
        for k=1:length(scale_small_probs)
            for l=1:length(data_nums)
                error=zeros(1,length(data_nums));
                error(l)=((MLE_q_approx_simulation(i,j,k,:))-(MLE_q_analytic));
                plot(data_nums,abs(error(l)))
            end
        end
    end
end
end
