function [P] = GraphMatcherPolynomial(G1,G2)
  a=size(G1);
b=size(G2);

if (a(1) == b(1))&&(a(2)==b(2));
    s = 0;
    N=a(1);
    d=[.65:.001:.88];
    while s<N
        s=s+1;
    %first constants
    fprintf('\n computing match in trial s = %i damping = %f \n',s,d);
    Onevec=ones(N,1);
    
    
    %first graph
    Dx1=sum(G1,2);
        for i = 1:N
            if Dx1(i)==0;
            Dx1(i)=1;
            end
        end
    
    thetax1=diag(Dx1.^-1);
    Wx1=(thetax1*G1)';
    X1(:,s)=((1-d(s))/N)*inv(eye(N)-d(s)*Wx1)*Onevec;
    
    %reverse the directed graph
     Dy1=sum(G1',2);
        for i = 1:N
            if Dy1(i)==0;
            Dy1(i)=1;
            end
        end
    
    thetay1=diag(Dy1.^-1);
    Wy1=(thetay1*G1')';
    Y1(:,s)=((1-d(s))/N)*inv(eye(N)-d(s)*Wy1)*Onevec;
    
    
    %graph two
    Dx2=sum(G2,2);
        for i = 1:N
            if Dx2(i)==0;
            Dx2(i)=1;
            end
        end
    thetax2=diag(Dx2.^-1);
    Wx2=(thetax2*G2)';
    X2(:,s)=((1-d(s))/N)*inv(eye(N)-d(s)*Wx2)*Onevec;
  %  disp([X1,X2]);
     Dy2=sum(G2',2);
    
        for i = 1:N
            if Dy2(i)==0;
            Dy2(i)=1;
            end
        end
    thetay2=diag(Dy2.^-1);
    Wy2=(thetay2*G2')';
    Y2(:,s)=((1-d(s))/N)*inv(eye(N)-d(s)*Wy2)*Onevec;
    
    %%comparing;
    P1xcat=horzcat(X1(:,s),eye(N));
    P1xcat=sortrows(P1xcat,1);
    P2xcat=horzcat(X2(:,s),eye(N));
    P2xcat=sortrows(P2xcat,1);
    P1ycat=horzcat(Y1(:,s),eye(N));
    P1ycat=sortrows(P1ycat,1);
    P2ycat=horzcat(Y2(:,s),eye(N));
    P2ycat=sortrows(P2ycat,1);
    
   
    
 
    P1realx=P1xcat(:,2:N+1);
    P2realx=P2xcat(:,2:N+1);
    Psupposedx=inv(P2realx)*P1realx;
    
  
    P1realy=P1ycat(:,2:N+1);
    P2realy=P2ycat(:,2:N+1);
    Psupposedy=inv(P2realy)*P1realy;
    

    
    if((G2==Psupposedx*G1*Psupposedx')&(G2==Psupposedy*G1*Psupposedy'))
        fprintf(' \n\n !!!!!!MATCH FOUND!!!!!\n\n\n\n Returning Permutation Matrix');
    P=Psupposedx;
    return
    end

       
    
       fprintf('no match found');
       P=zeros(N);
    end  
    fprintf('\n\n\n NO MATCH FOUND \n\n\n');
    
    
else
fprintf('functions have incorrect function size');
P=zeros(N);
end 
end 