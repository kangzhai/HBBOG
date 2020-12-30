function [Alpha_score,Cbest,Alpha_pos]=HBBOG(f,a,b,d,M,N)%混合差分扰动和趋优迁移的BBO算法
% Hybrid Biogeography Based Optimization with Grey wolf optimization(HBBOG)

Alpha_pos=zeros(1,d);
Alpha_score=inf; %change this to -inf for maximization problems
Beta_pos=zeros(1,d);
Beta_score=inf; %change this to -inf for maximization problems
Delta_pos=zeros(1,d);
Delta_score=inf; %change this to -inf for maximization problems
lambdaLower = 0.0; % lower bound for immigration probabilty per gene
lambdaUpper = 1; % upper bound for immigration probabilty per gene
Ib = 1; % max immigration rate for each island
Eb = 1; % max emigration rate, for each island
P = N; % max species count, for each island
dn=d+1;
X=rand(N,d+1);
aa=repmat(a,N,1);
bb=repmat(b,N,1);
X(:,1:d)=aa+(bb-aa).*X(:,1:d);
X(:,d+1)=feval(f, X(:,1:d));fit=X(:,d+1);
X = PopSort(X);
for i=1:N
    if fit(i)<Alpha_score
        Alpha_score=fit(i); % Update alpha
        Alpha_pos=X(i,1:d);
    end
    if fit(i)>Alpha_score && fit(i)<Beta_score
        Beta_score=fit(i); % Update beta
        Beta_pos=X(i,1:d);
    end
    if fit(i)>Alpha_score && fit(i)>Beta_score && fit(i)<Delta_score
        Delta_score=fit(i); % Update delta
        Delta_pos=X(i,1:d);
    end
end

SpeciesCount=zeros(1);
lambda=zeros(1);
mu=zeros(1);
for j = 1 : N
    if X(j,d+1) < inf
        SpeciesCount(j) = P - j;
    else
        SpeciesCount(j) = 0;
    end
    lambda(j) = Ib * (1 - SpeciesCount(j) / P);
    mu(j) = Eb * SpeciesCount(j) / P;    
end
lambdaMin = min(lambda);
lambdaMax = max(lambda);
lambdaScale = lambdaLower + (lambdaUpper - lambdaLower) * (lambda - lambdaMin) / (lambdaMax - lambdaMin);
Cbest=zeros(1,M);
Ped=0.5;a1=4;a2=2;
l=0;
while l<M    
    p=sqrt(l/M);
    Positions=X(:,1:d);
    as=1-l*((1-0)/M);
    if l<1*M/2
        for k = 1 : N
            alfa=rand^1;
            j=ceil(rand*d);            
            if rand < lambdaScale(k)
                % Pick a habitat from which to obtain a feature
                SelectIndex=ceil((k-1)*rand);                
                if rand>p
                    if rand<0.5
                        Positions(k,j) =X(SelectIndex,j);
                    else
                       cnum=ceil(rand*d);Positions(k,j) =X(SelectIndex,cnum); 
                    end
                else                    
                    Positions(k,j) =X(SelectIndex,j)+2*p*(0.5-rand)*(X(k,j)-X(SelectIndex,j));
                end
            else
                rnum=ceil(N*rand);while rnum==k,rnum=ceil(N*rand);end
                rnum1=ceil(N*rand);while rnum1==k ||rnum1==rnum,rnum1=ceil(N*rand);end                
                Positions(k,j) = X(k,j)+alfa*(X(1,j)-X(k,j)+X(rnum1,j)-X(rnum,j));                              
            end            
        end        
    else
        if rand<Ped
            Num=ceil(N*rand);           
            for i=1:N
                if i==Num                   
                    Positions(i,:)=a+(b-Alpha_pos);
                else
                    if rand<=0.5
                        k=ceil(rand*d);
                        X1=Alpha_pos(k)-a1*as*(2*rand-1)*abs(2*rand*Alpha_pos(k)-Positions(i,k));
                        X2=Beta_pos(k)-a1*as*(2*rand-1)*abs(2*rand*Beta_pos(k)-Positions(i,k));
                        X3=Delta_pos(k)-a1*as*(2*rand-1)*abs(2*rand*Delta_pos(k)-Positions(i,k));
                        Positions(i,k)=(X1+X2+X3)/3;% Equation (3.7)
                    else                        
                        X1=Alpha_pos-a2*as*(2*rand(1,d)-1).*abs(2*rand(1,d).*Alpha_pos-Positions(i,:));
                        X2=Beta_pos-a2*as*(2*rand(1,d)-1).*abs(2*rand(1,d).*Beta_pos-Positions(i,:));
                        X3=Delta_pos-a2*as*(2*rand(1,d)-1).*abs(2*rand(1,d).*Delta_pos-Positions(i,:));
                        Positions(i,:)=(X1+X2+X3)/3;
                    end
                    
                end
            end
        else            
            for k = 1 : N
                alfa=rand^1;
                for j = 1 : d
                    if rand < lambdaScale(k)
                        % Pick a habitat from which to obtain a feature
                        SelectIndex=ceil((k-1)*rand);                       
                        if rand>p
                            if rand<0.5
                            cnum=ceil(rand*d);Positions(k,j) =X(SelectIndex,cnum);
                            else
                                Positions(k,j) =X(SelectIndex,j);
                            end
                        else
                            Positions(k,j) =X(SelectIndex,j)+2*as*(0.5-rand)*(X(k,j)-X(SelectIndex,j));
                        end
                    else
                        rnum=ceil(N*rand);while rnum==k,rnum=ceil(N*rand);end
                        rnum1=ceil(N*rand);while rnum1==k ||rnum1==rnum,rnum1=ceil(N*rand);end                        
                        Positions(k,j) = X(k,j)+alfa*(X(1,j)-X(k,j)+X(rnum1,j)-X(rnum,j));                        
                    end
                end
            end            
        end        
    end   
    Positions=simplebounds(Positions,aa,bb);    
    fit=feval(f,Positions);    
    for i=1:N        
        fitness=fit(i);
        if X(i,d+1)>fitness
            X(i,d+1)=fitness;X(i,1:d)=Positions(i,:);
        end        
        % Update Alpha, Beta, and Delta
        if fitness<Alpha_score
            Alpha_score=fitness; % Update alpha
            Alpha_pos=Positions(i,:);
        end
        if fitness>Alpha_score && fitness<Beta_score
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
        end
        if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
        end        
    end
    X = PopSort(X);
    l=l+1;
    Cbest(l)=X(1,dn);
end
Alpha_score=X(1,dn);
Alpha_pos=X(1,1:d);

function s=simplebounds(s,Lb,Ub)
% Apply the lower bound
ns_tmp=s;
I=ns_tmp<Lb;
ns_tmp(I)=Lb(I);

% Apply the upper bounds
J=ns_tmp>Ub;
ns_tmp(J)=Ub(J);
% Update this new move
s=ns_tmp;