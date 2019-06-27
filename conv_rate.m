function[p,k_conv, SSE]=conv_rate(num_sims, e)
%e=avg_MLE_q_approx_simulation(:,i,j,k);
%for m=1:length(num_sims)


    y=log(1/e);
    x=log(num_sims);
    [p,c] = polyfit(x,y,1);
    %p=conv_fit(1);
    k_conv=exp(c);
    y_hat=polyval([p,k_conv],x);
    SSE=sum(y_hat-y).^2;
%end
end
   
    
    
    
    