function[c,p,sse] = convergence_rate_C(avg_MLE_error,num_sims)

F = @(x)error_function(avg_MLE_error,num_sims,x);
    x0=[0,1,1];
    
    [x,sse] = fmincon(F,x0,[],[],[],[],[0, -Inf],[Inf, Inf]);
    
    c=x(1);
    p=x(2);
    
end    

function [sse] = error_function(avg_MLE_error,num_sims,x)
    c=x(1);
    p=x(2);
    e_est = (1./(num_sims).^p) + c;
    sse = sum(((e_est) - (avg_MLE_error)).^2);
    
end