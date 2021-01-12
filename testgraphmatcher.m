
for N= 20:1000;
nodesize(N)=N;
A1 = randi([0 1], N,N);

b= randperm(N);

%not the same
A3 = randi([0 1], N,N);

P =zeros(N);
for i =1:N
P(i,b(i))=1;
end

P1=zeros(N);

A2=P*A1*P';
tic;
P1=GraphMatcherPolynomial(A1,A2);
T(N)=toc
if (P1==P);
    
fprintf('\nit works you are awesome!\n ');
    
end 
end 
plot(nodesize,T);