function[p,k, SSE]=conv_rate(num_sims, avg_MLE_q_approx_simulation(:,i,j,k))
e=avg_MLE_q_approx_simulation(:,i,j,k);
%for m=1:length(num_sims)


    y=log(1/e);
    x=log(num_sims(m));
    conv_fit = polyfit(x,y,1);
    p=conv_fit(1);
    k=e^conv_fit(2);
    y_hat=polyval(conv_fit,x);
    SSE=sum(y_hat-y).^2;
end
   
    
    
    
    