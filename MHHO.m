clc
close all
tic
%% Map grid
G=[0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0; 
   0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0; 
   0 1 1 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 1 1 1 0 0 0 0 0 1 1 1 1 0 0; 
   0 0 0 0 0 0 1 1 1 0 0 0 0 0 1 1 1 1 0 0; 
   0 1 1 0 0 0 1 1 1 0 0 0 0 0 1 1 1 0 0 0; 
   0 1 1 0 0 0 1 1 1 0 0 0 0 0 1 1 1 0 0 0;
   0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0; 
   0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0; 
   0 1 1 0 0 0 1 1 1 0 0 1 1 0 0 0 0 0 0 0; 
   0 1 1 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 1 1 1 1 0; 
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0; 
   1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0; 
   1 1 0 0 0 0 1 1 0 0 0 1 0 0 0 0 0 0 0 0; 
   0 0 0 0 0 0 1 1 0 1 1 1 0 0 0 0 0 1 1 0; 
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0; 
   0 0 1 1 0 0 0 0 0 0 1 1 0 0 1 0 0 0 0 0; 
   0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0;];

for i=1:20/2
    for j=1:20
        m=G(i,j);
        n=G(21-i,j);
        G(i,j)=n;
        G(21-i,j)=m;
    end
end
%% 
S = [1 1];   
E = [20 20];  
G0 = G;
G = G0(S(1):E(1),S(2):E(2)); 
[Xmax,dimensions] = size(G);        
dimensions = dimensions - 2;             

%% Parameter Setting
max_gen = 100;        % Maximum iteration number
num_polution = 20;    % Individuals of the population
X_min = 1;  

%% Initialization
X = zeros(num_polution,dimensions);
for i = 1:num_polution
    for j = 1:dimensions
       column = G(:,j+1);      
       id = find(column == 0); 
       X(i,j) =  id(randi(length(id))); 
       id = [];
    end 
    fit( i ) = fitness(X( i, : ),G);
end
fit_person_best = fit;  
person_best = X;     
[fit_global_best, best_person_Index] = min( fit );        
global_best = X(best_person_Index, : );   
[fit_max,B]=max(fit);
worse_person= X(B,:);  
%%
for gene = 1:max_gen
   gene
   [ans1,sort_Index] = sort(fit);   %适应值度从小到大排序
    [fit_max,B] = max(fit);
    worse_person = X(B,:);           %找出适应度最差的个体
    [~,Index] = sort(fit_person_best);
        
     for i=1:num_polution
      E1=2*(1-(gene/max_gen));  
      E0=2*rand()-1; %-1<E0<1
      Escaping_Energy=E1*(E0);  % escaping energy of rabbit   
         if abs(Escaping_Energy)>=1
            %% Exploration:
            % Harris' hawks perch randomly based on 2 strategy:
            q=rand();
            rand_Hawk_index = randi(length(num_polution));
            X_rand = person_best(Index(rand_Hawk_index),:);
            if q<0.5
                % perch based on other family members
                X(Index(i),:)= X_rand-rand()*abs(X_rand-2*rand()*person_best(Index(i),:));
                X(Index(i),:) = Bounds(X(Index(i),:), X_min,Xmax);                     
                fit(Index(i)) = fitness(X(Index(i),:),G);
            elseif q>=0.5
                % perch on a random tall tree (random site inside group's home range)
                 X(Index(i),:)=(global_best-mean(person_best(Index(i),:)))-rand()*((Xmax-X_min)*rand+X_min);
                 X(Index(i),:) = Bounds(X(Index(i),:), X_min,Xmax);                     
                 fit(Index(i)) = fitness(X(Index(i),:),G);
            end
          elseif abs(Escaping_Energy)<1
            %% Exploitation:
            % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit
            
            %% phase 1: surprise pounce (seven kills)
            % surprise pounce (seven kills): multiple, short rapid dives by different hawks
            
            r=rand(); % probablity of each event
            
            if r>=0.5 && abs(Escaping_Energy)<0.5 % Hard besiege 
               X(Index(i),:)=global_best-Escaping_Energy*abs(global_best-person_best(Index(i),:));
               X(Index(i),:) = Bounds(X(Index(i),:), X_min,Xmax);                     
                fit(Index(i)) = fitness(X(Index(i),:),G);
            end
            
            if r>=0.5 && abs(Escaping_Energy)>=0.5  % Soft besiege
                Jump_strength=2*(1-rand()); % random jump strength of the rabbit
                X(Index(i),:)=(global_best-person_best(Index(i),:))-Escaping_Energy*abs(Jump_strength*global_best-person_best(Index(i),:));
                X(Index(i),:) = Bounds(X(Index(i),:), X_min,Xmax);                     
                fit(Index(i)) = fitness(X(Index(i),:),G);
            end  
         %% phase 2: performing team rapid dives (leapfrog movements)
            if r<0.5 && abs(Escaping_Energy)>=0.5 % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                w1=2*exp(-(8*gene/max_gen)^2);
                Jump_strength=2*(1-rand());
                X1(Index(i),:)=w1*global_best-Escaping_Energy*abs(Jump_strength*global_best-person_best(Index(i),:));
                X1(Index(i),:) = Bounds(X1(Index(i),:), X_min,Xmax);                      
               fit1(Index(i)) = fitness(X1(Index(i),:),G);
                if  fit1(Index(i))< fit(Index(i)) % improved move
                     X(Index(i),:)= X1(Index(i),:);
                else % hawks perform levy-based short rapid dives around the rabbit
                     X2(Index(i),:)=w1*global_best-Escaping_Energy*abs(Jump_strength*global_best-person_best(Index(i),:))+rand(1,dimensions).*Levy(dimensions);
                     X2(Index(i),:) = Bounds(X2(Index(i),:), X_min,Xmax);                     
                     fit2(Index(i)) = fitness(X2(Index(i),:),G);
                    if (fit2(Index(i))<fit2(Index(i))) % improved move?
                        X(Index(i),:)= X2(Index(i),:);
                    end
                end
            end
            
            if r<0.5 && abs(Escaping_Energy)<0.5 % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                % hawks try to decrease their average location with the rabbit
                w1=2*exp(-(8*gene/max_gen)^2);
                Jump_strength=2*(1-rand()); 
                X1(Index(i),:)=w1*global_best-Escaping_Energy*abs(Jump_strength*global_best-mean(person_best(Index(i),:)));
                X1(Index(i),:) = Bounds(X1(Index(i),:), X_min,Xmax);         
                fit1(Index(i)) = fitness(X1(Index(i),:),G);
                
                if fit1(Index(i))<fit(Index(i)) % improved move?
                     X(Index(i),:) = X1(Index(i),:) ;
                else % Perform levy-based short rapid dives around the rabbit 
                     X2(Index(i),:)=w1*global_best-Escaping_Energy*abs(Jump_strength*global_best-mean(person_best(Index(i),:)))+rand(1,dimensions).*Levy(dimensions);
                     X2(Index(i),:) = Bounds(X2(Index(i),:), X_min,Xmax);     
                     fit2(Index(i)) = fitness(X2(Index(i),:),G);
                    if ( fit2(Index(i))< fit(Index(i))) % improved move?
                        X(Index(i),:)=X2(Index(i),:);
                    end
                end
            end             
         end 
     end
    
     for i=1:num_polution
         X(i,:) = Bounds(X(i,:),X_min,Xmax);
     end
     for i = 1:num_polution
        if (fit(i)<fit_person_best(i))
            fit_person_best(i) = fit(i);
            person_best(i,:) = X(i,:);
        end
        if(fit_person_best(i) < fit_global_best)
            fit_global_best = fit_person_best(i);
            global_best = person_best(i,:);
        end
    end
    global_best = LocalUpdateSearch(global_best,Xmax,G);%%%%%Local update search
    fit_global_best = fitness(global_best,G);
    final_goal(gene,1)=fit_global_best;
end
toc
%% result 
global_best = round(global_best);
fit_global_best
figure(1)
plot(final_goal,'b-','LineWidth',1.2);
xlabel('Iteration number');
ylabel('Optimum path');
route = [S(1) global_best E(1)];
path=GenerateRoute(route,G);
path=GenerateSmoothPath(path,G);
path=GenerateSmoothPath(path,G);
figure(2)
for i=1:20/2
    for j=1:20
        m=G(i,j);
        n=G(21-i,j);
        G(i,j)=n;  
        G(21-i,j)=m;
    end
end  
n=20;
for i=1:20
    for j=1:20
        if G(i,j)==1 
            x1=j-1;y1=n-i; 
            x2=j;y2=n-i; 
            x3=j;y3=n-i+1; 
            x4=j-1;y4=n-i+1; 
            fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.14,0.87,0.87]); 
            hold on 
        else 
            x1=j-1;y1=n-i; 
            x2=j;y2=n-i; 
            x3=j;y3=n-i+1; 
            x4=j-1;y4=n-i+1; 
            fill([x1,x2,x3,x4],[y1,y2,y3,y4],[1,1,1]); 
            hold on 
        end 
    end 
end 
hold on
xlabel('Environment 1')
drawPath(path,G)
title('MHHO')

function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end

