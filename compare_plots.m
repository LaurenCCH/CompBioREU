function[]=compare_plots(MLE_q_approx_simulation,MLE_q_analytic)
%This function loops over ordered triples (data_nums,bw,scale_small_probs)and
%plots num_sims verusus abs(MLE_q_approx_simulation-MLE_q_analytic).