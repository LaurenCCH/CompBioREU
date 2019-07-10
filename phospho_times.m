function [TT]=phospho_times(q,h,n)
%This code tracks the time it takes proteins to be phosphorylated.
%n is the number of protiens that are tracked
%we assume the probability a protein is phosphorylated over a time step of
%length h is poisson distributed with m=1 (the number of events is 1).
hold off

%memoizedFactorial = memoize(@factorial);

%p=@(q,h,m)((exp(-q*h).*(q*h).^m)/memoizedFactorial(m));
%p=@(q,h,m)(q*h);

%p = q*h;
t=zeros(1,n);
%m is the number of phosphorylations.  It must be set to 1.
%m=1;

n=int64(n);
chunk_size=int64(100);

num_chunks=ceil(n/chunk_size);

size_last_chunk=mod(n,chunk_size);

num_full_chunks=floor(n/chunk_size);

chunk_sizes=zeros(1,num_chunks);
chunk_sizes(1:num_full_chunks)=chunk_size;

if size_last_chunk~=0
    chunk_sizes(length(chunk_sizes))=size_last_chunk;
end
t=zeros(num_chunks,chunk_size);
T=cell(num_chunks);
TT=zeros(1,n);
parfor chunk_index=1:num_chunks
    t=zeros(1,chunk_sizes(chunk_index));
    for s=1 : chunk_sizes(chunk_index)
    %r=rand;
    i=0;
    %while rand > p(q,h,m)
    while rand > q*h
        i=i+1;
        %r=rand;
    end
    t(s) = i*h;
    end
    T{chunk_index}=t;
end

for chunk_index=1:num_chunks
    TT(chunk_size*(chunk_index-1)+1:chunk_size*(chunk_index-1)+chunk_sizes(chunk_index))=T{chunk_index};
end


histogram(TT,'Normalization','pdf')
hold on
tau=0:.1:20;
plot(tau,q.*exp(-q*tau))
hold off
end

%Need a liklihood of the data, find q that maximizes the liklihood
%try analytic, or find liklihood function and evaluate numericaly with
%matlab
  