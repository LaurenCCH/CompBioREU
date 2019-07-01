function[p,k, SSE]=convergence_rate(num_sims,avg_MLE_error)
%we expect avg_MLE_error to converge to 0 like (k/num_sims^p) as num_sims
%grows.
%reparameterizing into x and y for a linear relationship.
    
y=log(1./avg_MLE_error)';
x=log(num_sims);
params = polyfit(x,y,1);
p=params(1);
c=params(2);
k=exp(c);
y_hat=polyval([p,c],x);
SSE=sum((y_hat-y).^2);
end