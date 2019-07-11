function [c,k,p] = convergence_rate_C(avg_MLE_error,num_sims)

F = @(c,k,p)error_function(avg_MLE_error,num_sims,c,k,p);
x0=[0,1,1];
    [x,sse] = fminsearch(F,x0);
    x(1) = c;
    x(2) = k;
    x(3) = p;
end    

function [sse] = error_function(avg_MLE_error,num_sims,c,k,p)

    e_est = (k/(num_sims).^p) + c;
    sse = sum(((e_est) - (avg_MLE_error)).^2);
    
end