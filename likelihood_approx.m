%This code take a time and a q as inputs and uses the generated data to
%calculate the probability of observing a single event within a range of
%times

function [log_liklihood_t_values_1, list_probs] = likelihood_approx(t,q0,num_sims,bw,degenerate_probability)


%Here we set the number of simulations as the current iteration of num_sims
n=num_sims;

% recalling phospho_times to calculate the times for q0
[t0] = phospho_times(q0,.01,n);

%This code creates a histogram and then takes the bin counts and turns
%them into probabilities
histogram(t0,'BinWidth', bw,'Normalization','pdf')

%This code take the bin counts and turns them into probabilities
[N,edges] = histcounts(t0,'BinWidth', bw,'Normalization', 'probability');

%This divides the probabilities by the bin-width to calculate the height of
%the bin
N_1 = N./bw;


%This puts the vector of times into numerical order
t_1=sort(t);

%This initializes a vector, list_probs, with zeros
list_probs=zeros(1,length(t));

%This loop iterates through the vector of sorted times, t_1, and 
%sorts each time into the proper bin and thus assigns it the corresponding
%probability
for d=1:length(t)
    %If the time is less than the left bound of the histogram
    % or greater than the right bound, it assigns a probability of zero
    if t_1(d)<edges(1) || t_1(d)>edges(length(edges))
        list_probs(d)=0;
    
    else
    
        %If the time is on the rightmost edge of the histogram, it is
        %assigned the probability corresponding to the last bin
        if t_1(d)==edges(length(edges))
            index_t = floor(t_1(d)/bw);
        
        %Otherwise, the time value is divided by the bin width, rounded down,
        %and one is added to it to find which index of N_1 corresponds to the
        %time.
        else
            index_t = floor(t_1(d)/bw)+1;
        end
        
    %Once the time has been matched to the proper probability, that 
    %probability is put into the list_probs vector at the index that 
    %matches the time index
    list_probs(d) = N_1(index_t);
    end
end
       
%This loop assigns a degenarte probability to every probability of zero to
%prevent the log likelihood from equaling -Inf
for m=1:length(list_probs)
    if list_probs(m)==0
        list_probs(m)=degenerate_probability;
    end
end

%This takes the log of every probability and sums them to calculate the log
%likelihood of t.
log_liklihood_t_values_1 = sum(log(list_probs));


end


	
	
