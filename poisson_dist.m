function [x]=euler_method(m,q,t,n)
%This code tracks the time it takes proteins to be phosphorylated.
%n is the number of protiens that are tracked
%we assume the probability a protein is phosphorylated over a time step of
%length t is poisson distributed with m=1 (the number of events is 1).
hold off

p=@(q,t,m)((exp(-q*t)*(q*t)^m)/factorial(m));
x=zeros(1,n);
for s=1 : n
    
    %r=rand;
    i=0;
    while rand > p(q,t,m)
        i=i+1;
        %r=rand;
    end
   
    x(s) = i*t;
    
   
end
histogram(x,'Normalization','pdf')
hold on
tau=0:.1:20;
plot(tau,exp(-q*tau).*(q))
hold off
end

%Need a liklihood of the data, find q that maximizes the liklihood
%try analytic, or find liklihood function and evaluate numericaly with
%matlab
  