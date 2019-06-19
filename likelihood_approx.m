%This code take a time and a q as inputs and uses the generated data to
%calculate the probability of observing a single event within a range of
%times

%This code take a time and a q as inputs and uses the generated data to
%calculate the probability of observing a single event within a range of
%times
function [log_liklihood_t_values_1, list_probs] = likelihood_approx(t,q0)
% recalling phospho_times to calculate the times for q0 
bw=.1;
n=10000;
[t0] = phospho_times(q0,.01,n);
%This code creates a histogram and then takes the bin counts and turns 
%them into probabilities
histogram(t0,'BinWidth', bw,'Normalization','pdf');
[N,edges] = histcounts(t0,'BinWidth', bw,'Normalization', 'probability');
N_1 = N./bw;
sum(N)
saveas(gcf,fprintf('likelihood_plot_approx%f', round(n)))
%This code take the bin counts and turns them into probabilities

%This initializes A1 with zeros
A1=zeros(1,((round(edges(end))*3)+1));
r=0;
b=1;
%This fills A1 with the time boundaries and the associated probabilities
% for a=1:(edges(end))
% 
% 
%    
%     
%     A1(a)=r;
%     
%     r=r+bw;
%     b=b+1;
%     if b>length(N)
%     break
%     end
%     
% end

t_1=sort(t);
N_1
list_probs=zeros(1,length(t));

for d=1:length(t)
    index_t = ceil(t_1(d)/bw);
    list_probs(d) = N_1(index_t);
end
    
    
    
   
    
%     while (t_1(d)>=A1(c+1) && c+1<=length(N))
%        % if c>length(A1)
%         %break
%         %end
%         c=c+1;
%         if c==length(A1)
%             break
%         end
%     end
%     list_probs(d)=N(c);
        
    


        
        
    
%This code formats the data into an easily readable format
%formatSpec = 'The probability for %2.4f <t< %2.4f is %1.5f \n' ;
%fprintf(formatSpec,A1)  
    

%log_liklihood_t_values = sum(log(list_probs));
for m=1:length(list_probs)
    if list_probs(m)==0
        list_probs(m)=(1/(n*10));
    else
        list_probs(m)=list_probs(m);
    end
end
%list_probs
log_liklihood_t_values_1 = sum(log(list_probs));


end
