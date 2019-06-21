function [t]=phospho_times(q,h,n)
%This code tracks the time it takes proteins to be phosphorylated.
%n is the number of protiens that are tracked
%we assume the probability a protein is phosphorylated over a time step of
%length h is poisson distributed with m=1 (the number of events is 1).
hold off

%memoizedFactorial = memoize(@factorial);

%p=@(q,h,m)((exp(-q*h).*(q*h).^m)/memoizedFactorial(m));
%p=@(q,h,m)(q*h);


t=zeros(1,n);
%m is the number of phosphorylations.  It must be set to 1.
m=1;
for s=1 : n
    
    %r=rand;
    i=0;
    %while rand > p(q,h,m)
    while rand > q*h
        i=i+1;
        %r=rand;
    end
   
    t(s) = i*h;
    
   
end
histogram(t,'Normalization','pdf')
hold on
tau=0:.1:20;
plot(tau,q.*exp(-q*tau))
hold off
end

%Need a liklihood of the data, find q that maximizes the liklihood
%try analytic, or find liklihood function and evaluate numericaly with
%matlab
  