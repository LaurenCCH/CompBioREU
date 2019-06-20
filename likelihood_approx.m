%This code take a time and a q as inputs and uses the generated data to
%calculate the probability of observing a single event within a range of
%times

%This code take a time and a q as inputs and uses the generated data to
%calculate the probability of observing a single event within a range of
%times
function [log_liklihood_t_values_1, list_probs] = likelihood_approx(t,q0)
% recalling phospho_times to calculate the times for q0 
bw=.01;
n=1000;
[t0] = phospho_times(q0,.01,n);
%This code creates a histogram and then takes the bin counts and turns 
%them into probabilities
histogram(t0,'BinWidth', bw,'Normalization','pdf');
[N,edges] = histcounts(t0,'BinWidth', bw,'Normalization', 'probability');
N_1 = N./bw;

%saveas(gcf,fprintf('likelihood_plot_approx%f', round(n)))
%This code take the bin counts and turns them into probabilities

t_1=sort(t);

list_probs=zeros(1,length(t));

for d=1:length(t)
    if t_1(d)<edges(1) || t_1(d)>edges(length(edges))
        list_probs(d)=0;
    else
    index_t = floor(t_1(d)/bw)+1;
    list_probs(d) = N_1(index_t);
    end
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
    end
end
%list_probs
log_liklihood_t_values_1 = sum(log(list_probs));


end
