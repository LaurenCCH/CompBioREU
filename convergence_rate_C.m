function[c,k,p,sse] = convergence_rate_C(avg_MLE_error,num_sims)

F = @(c,k,p)error_function(avg_MLE_error,num_sims,c,k,p);
    x0=[0,1,1];
    [x,sse] = fminsearch(F,x0);
    
    c=x(1);
    k=x(2);
    p=x(3);
    
end    

function [sse] = error_function(avg_MLE_error,num_sims,c,k,p)

    e_est = (k/(num_sims).^p) + c;
    sse = sum(((e_est) - (avg_MLE_error)).^2);
    
end