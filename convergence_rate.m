function[p,k, SSE]=convergence_rate(num_sims,avg_MLE_error)
%we expect avg_MLE_error to converge to 0 like (k/num_sims^p) as num_sims
%grows.
%reparameterizing into x and y for a linear relationship.
    
y=log(avg_MLE_error);
x=log(num_sims);
params = polyfit(x,y,1);
p=-params(1);
c=params(2);
k=exp(c);
y_hat=polyval([params(1),params(2)],x);
avg_MLE_error_hat=exp(y_hat);
%SSE=sum((y_hat-y).^2);
SSE=sum((avg_MLE_error_hat-avg_MLE_error).^2);

%plot(x,y,'*',x,y_hat,'o')
end