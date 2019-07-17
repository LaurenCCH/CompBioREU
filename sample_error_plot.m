function[]=sample_error_plot(sample_error_hist,data_nums,bw,scale_small_probs, num_sims)
histogram(sample_error_hist,20)
annotation('textbox','String',sprintf("data_nums=%.2f\nbw=%.2f\nscale_small_probs=%.2f\nnum_sims=%.2f", data_nums,bw,scale_small_probs, num_sims),'FitBoxToText','on', 'Interpreter', 'none', 'Position', [.5, .7, .4, .2]);
xlabel('Sample error','FontSize', 16)
ylabel('Frequency', 'FontSize', 16)

end