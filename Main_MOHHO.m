

clc;
clear;
close all;

%% Problem Definition

CostFunction=@(x) CostFunctio(x);      % Cost Function

nVar=5;             % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=0;          % Lower Bound of Variables
VarMax=1;          % Upper Bound of Variables


%% MOPSO Parameters

MaxIt=200;           % Maximum Number of Iterations

nPop=200;            % Population Size

nRep=50;            % Repository Size

w=0.5;              % Inertia Weight
wdamp=0.99;         % Intertia Weight Damping Rate
c1=1;               % Personal Learning Coefficient
c2=2;               % Global Learning Coefficient

nGrid=7;            % Number of Grids per Dimension
alpha=0.1;          % Inflation Rate

beta=2;             % Leader Selection Pressure
gamma=2;            % Deletion Selection Pressure

mu=0.1;             % Mutation Rate

%% Initialization

empty_Rabbit.Location=[];
empty_Rabbit.Cost=[];
empty_Rabbit.Sol=[];
empty_Rabbit.IsDominated=[];
empty_Rabbit.GridIndex=[];
empty_Rabbit.GridSubIndex=[];

Rabbits=repmat(empty_Rabbit,nPop,1);
X = zeros(nPop, nVar);
Rabbit_Location=zeros(VarSize);
Rabbit_Energy=inf;

for i=1:nPop
    
    Rabbits(i).Location = rand(VarSize).*(VarMax-VarMin)+VarMin; 
    X(i,:) = rand(VarSize).*(VarMax-VarMin)+VarMin; 
    [Rabbits(i).Cost, Rabbits(i).Sol] = CostFunction(Rabbits(i).Location);
    
end

% Determine Domination
Rabbits=DetermineDomination(Rabbits);

rep=Rabbits(~[Rabbits.IsDominated]);

Grid=CreateGrid(rep,nGrid,alpha);

for i=1:numel(rep)
    rep(i)=FindGridIndex(rep(i),Grid);
end


%% MOPSO Main Loop

for it=1:MaxIt
    E1=2*(1-(it/MaxIt)); % factor to show the decreaing energy of rabbit    
    for i=1:nPop
        
        leader=SelectLeader(rep,beta);
        
        
        E0=2*rand()-1; %-1<E0<1
        Escaping_Energy=E1*(E0);  % escaping energy of rabbit
        
        if abs(Escaping_Energy)>=1
            %% Exploration:
            % Harris' hawks perch randomly based on 2 strategy:
            
            q=rand();
            rand_Hawk_index = floor(nPop*rand()+1);
            X_rand = Rabbits(rand_Hawk_index);
            if q<0.5
                % perch based on other family members
                Rabbits(i).Location=X_rand.Location-rand()*abs(X_rand.Location-2*rand()*Rabbits(i).Location);
                X(i,:)=X_rand.Location-rand()*abs(X_rand.Location-2*rand()*Rabbits(i).Location);
            elseif q>=0.5
                % perch on a random tall tree (random site inside group's home range)
                Rabbits(i).Location=(leader.Location-mean(X))-rand()*((VarMax-VarMin)*rand+VarMin);
                X(i,:)=(leader.Location-mean(X))-rand()*((VarMax-VarMin)*rand+VarMin);
            end
            
        elseif abs(Escaping_Energy)<1
            %% Exploitation:
            % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit
            
            %% phase 1: surprise pounce (seven kills)
            % surprise pounce (seven kills): multiple, short rapid dives by different hawks
            
            r=rand(); % probablity of each event
            
            if r>=0.5 && abs(Escaping_Energy)<0.5 % Hard besiege
                Rabbits(i).Location=(leader.Location)-Escaping_Energy*abs(leader.Location-Rabbits(i).Location);
                X(i,:)=(leader.Location)-Escaping_Energy*abs(leader.Location-X(i,:));
            end
            
            if r>=0.5 && abs(Escaping_Energy)>=0.5  % Soft besiege
                Jump_strength=2*(1-rand()); % random jump strength of the rabbit
                X(i,:)=(leader.Location-X(i,:))-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
                Rabbits(i).Location=(leader.Location-Rabbits(i).Location)-Escaping_Energy*abs(Jump_strength*Rabbit_Location-Rabbits(i).Location);
            end
            
            %% phase 2: performing team rapid dives (leapfrog movements)
            if r<0.5 && abs(Escaping_Energy)>=0.5 % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                
                Jump_strength=2*(1-rand());
                X1.Location=leader.Location-Escaping_Energy*abs(Jump_strength*leader.Location-X(i,:));
                [X1.Cost, X1.Sol] = CostFunction(X1.Location);
                if Dominates(X1,Rabbits(i))
                    Rabbits(i).Location=X1.Location;
                    Rabbits(i).Cost=X1.Cost;
                    Rabbits(i).Cost=X1.Cost;

                elseif Dominates(Rabbits(i),X1)
                    X2.Location=leader.Location-Escaping_Energy*abs(Jump_strength*leader.Location-X(i,:))+rand(1,nVar).*Levy(nVar);                    
                    [X2.Cost, X2.Sol] = CostFunction(X2.Location);
                    if Dominates(X2,Rabbits(i))
                        Rabbits(i).Location=X2.Location;
                        Rabbits(i).Cost=X2.Cost;
                        Rabbits(i).Cost=X2.Cost;
                    end
                else
                    if rand<0.5
                        Rabbits(i).Location=X1.Location;
                        Rabbits(i).Cost=X1.Cost;
                        Rabbits(i).Cost=X1.Cost;
                        Rabbits(i).Sol=X1.Sol;
                    end
                end                                
            end
            
            if r<0.5 && abs(Escaping_Energy)<0.5 % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                % hawks try to decrease their average location with the rabbit
                Jump_strength=2*(1-rand());
                X1.Location=leader.Location-Escaping_Energy*abs(Jump_strength*leader.Location-mean(X));
                [X1.Cost, X1.Sol] = CostFunction(X1.Location);
                if Dominates(X1,Rabbits(i))
                    Rabbits(i).Location=X1.Location;
                    Rabbits(i).Cost=X1.Cost;
                    Rabbits(i).Cost=X1.Cost;

                elseif Dominates(Rabbits(i),X1)
                    X2.Location=leader.Location-Escaping_Energy*abs(Jump_strength*leader.Location-mean(X))+rand(1,nVar).*Levy(nVar);                    
                    [X2.Cost, X2.Sol] = CostFunction(X2.Location);
                    if Dominates(X2,Rabbits(i))
                        Rabbits(i).Location=X2.Location;
                        Rabbits(i).Cost=X2.Cost;
                        Rabbits(i).Cost=X2.Cost;
                    end
                else
                    if rand<0.5
                        Rabbits(i).Location=X1.Location;
                        Rabbits(i).Cost=X1.Cost;
                        Rabbits(i).Cost=X1.Cost;
                        Rabbits(i).Sol=X1.Sol;
                    end
                end                               
            end
        end
        
        Rabbits(i).Location = max(Rabbits(i).Location, VarMin);
        Rabbits(i).Location = min(Rabbits(i).Location, VarMax);
        
        [Rabbits(i).Cost, Rabbits(i).Sol] = CostFunction(Rabbits(i).Location);
        
        % Apply Mutation
        pm=(1-(it-1)/(MaxIt-1))^(1/mu);
        if rand<pm
            NewSol.Location=Mutate(Rabbits(i).Location,pm,VarMin,VarMax);
            [NewSol.Cost, NewSol.Sol]=CostFunction(NewSol.Location);
            if Dominates(NewSol,Rabbits(i))
                Rabbits(i).Location=NewSol.Location;
                Rabbits(i).Cost=NewSol.Cost;
                Rabbits(i).Sol=NewSol.Sol;

            elseif Dominates(Rabbits(i),NewSol)
                % Do Nothing

            else
                if rand<0.5
                    Rabbits(i).Location=NewSol.Location;
                    Rabbits(i).Cost=NewSol.Cost;
                    Rabbits(i).Sol=NewSol.Sol;
                end
            end
        end
        
    end
    
    % Add Non-Dominated Particles to REPOSITORY
    rep=[rep
         Rabbits(~[Rabbits.IsDominated])]; %#ok
    
    % Determine Domination of New Resository Members
    rep=DetermineDomination(rep);
    
    % Keep only Non-Dminated Memebrs in the Repository
    rep=rep(~[rep.IsDominated]);
    
    % Update Grid
    Grid=CreateGrid(rep,nGrid,alpha);

    % Update Grid Indices
    for i=1:numel(rep)
        rep(i)=FindGridIndex(rep(i),Grid);
    end
    
    % Check if Repository is Full
    if numel(rep)>nRep
        
        Extra=numel(rep)-nRep;
        for e=1:Extra
            rep=DeleteOneRepMemebr(rep,gamma);
        end
        
    end
    
    % Plot Costs
    figure(1);
    PlotCosts(Rabbits,rep);
    pause(0.01);
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);
    
    % Damping Inertia Weight
    w=w*wdamp;
    
end

%% Resluts

solutions = [];
costs = [];
for i=1:numel(rep)
solutions = cat(1, solutions, rep.Location(i));
costs = cat(1, costs, rep.Cost(i));
Sols = cat(1, Sols, rep.Sol(i));
end

