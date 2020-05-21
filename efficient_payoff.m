%set up
N = 100; % the size of the social network
A = zeros(N,N); % initial adjacency matrix

doTest = true;

alpha0 = 3; % the direct benefit
cost = 1; % the cost of a link

%compute the initial payoff of each individual in the network
[payoff, A2] = structHolePayoff(A,alpha0,cost); 

noswaps=700; % make 1000 edge swaps
t = 1;
 while (t<=noswaps)
     [payoff,A,A2,time_efficient, time_inefficient, edges] = ...
         efficientPayoff(A,alpha0,cost,t);
 
     t=t+1;
 end
     A = double((A+A')~=0); % make the network undirected
     % plot the network
     G = graph(A,'omitselfloops');
     p = plot(G,'Layout','force','EdgeAlpha',0.005,'NodeColor','r');
     title('Generated network graph')
     
     figure(2)
     % plot the time of efficient and inefficient time
     plot(edges,time_inefficient)
     title('Efficient vs Inefficient method')
     hold on
     plot(edges,time_efficient)
     xlabel('edges')
     ylabel('time')
     legend('inefficient','efficient')
     
  [payoff1,t1] = efficientPayoff2(A);
  [payoff2,t2] = inefficientPayoff2(A);
  saved_time = (t2-t1)/t2;
  fprintf('The efficient method saved %0.3f%% \n', saved_time*100);
  
%%
function [payoff, A2] = structHolePayoff(A,alpha0,cost)
% A = adjacency matix
% alpha0 = benefit for a direct link
% cost = cost matrix, such that cost(u,v) is cost of forming link between
% nodes u and v
N=size(A,1);
if (nargin<2)
    % default value for direct benefit is 1.0
    alpha0=1.0;
end

if (nargin<3)
    % default value for cost is 1.0
    cost = ones(N,N)-sparse(1:N,1:N,ones(N,1));
end

% compute a matrix A2, such that A2(u,v)=number of length 2 
% paths between nodes u and v.  Set A2(u,v)=0, if A(u,v)=1, i.e.
% if the nodes are already directly connected,don't count their length 2
% paths

A2=((A*A).*(A==0));
A2 = A2 - sparse(1:N,1:N,diag(A2));
D = sum(A,2);

% compute Kleinberg's Harmonic benefit value for each node
% in the network
beta = harmonicBeta(A,A2,alpha0);

% add in the benefit alpha0 associated with all the direct links
% subtract away the cost assocatied with forming a link

payoff = alpha0*D-sum(cost.*A,2) + beta;

% add all the direct-link benefit alpha0 and the intermediary benefits
% subtract the cost associated with forming a link
% D is the number of direct links
% sum(matrix,2) adds all the numbers in the same row

end

function beta = harmonicBeta(A,A2,alpha0)
% A = adjacency matix
% A2 = matrix of length 2 paths.
% alpha0 = benefit associated with a direct link

N = size(A,1);
beta = zeros(N,1);

% run through each node v, and compute the benefit it
% obtains for being on a length 2 path between two nodes.

% To find the two nodes whose path v is on, get all the 
% neighbours directly connected to v, and read from A2 the
% total number of length two paths that connect two neighbours
for v=1:N
    mask = A(v,:)~=0; % mask to extract neighbours
    Av = A2(mask,mask); % sub-matrix of length 2 paths connecting neighbours
    
    % h is the benefit that node v obtains, which is the sum of the 
    % benefits over all pairs of neighbours connected by length-2 paths
    h = sum(sum(alpha0./Av(Av~=0)));
    beta(v) = h ;
end
end

%%
function [payoff,A,M,time_efficient, time_inefficient, edges] ...
    = efficientPayoff(A, alpha0, cost,t)
N = size(A,1);

%% STEP 1 
% Calculate the Payoff in the direct way from the full matrix

% Firstly calculate A2
A2=((A*A).*(A==0));
A2 = A2 - sparse(1:N,1:N,diag(A2));

% Secondly calculate the payoff
G = (1-A).*(A2~=0)./(A2+(A2==0));
H = A*G';
pold = sum(A'.*H,2);

% pold = zeros(N,1);
% for w=1:N
%     G = A(:,w).*A(w,:).*(1-A).*(A2~=0)./(A2+(A2==0));
%     pold(w) = sum(sum(G));
% end  
A2old = A2;

%% STEP 2
% Choose some pair (u0,v0) into which to put an edge
found=false;
while (~found)
    u0 = ceil(rand*N);
    p = randperm(N);
    for v=1:N
        v0 = p(v);
        if (v0==u0)
            continue;
        end
        if (A(u0,v0)==0)
            % Add the edge
            fprintf('Adding an edge at position (%d,%d)\n', u0,v0);
            A(u0,v0)=1;
            found=true;
            break;
        else
            % Remove the edge
            fprintf('Removing an edge at position (%d,%d)\n', u0,v0);
            A(u0,v0)=0;
            found=true;
            break;
        end
        
    end
end
%% STEP 3
% Get the new payoff in the direct way %%%%%%%
% Firstly calculate A2
tic
A2new=((A*A).*(A==0));
A2new = A2new - sparse(1:N,1:N,diag(A2new));

% Secondly calculate the payoff
G = (1-A).*(A2new~=0)./(A2new+(A2new==0));
H = A*G';
pnew = sum(A'.*H,2); 
time_inefficient(t) = toc;
edges(t) = sum(sum(A));
fprintf('time taken %d\n', time_inefficient(t));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 4
% Now update A2, only updating the rows and columns that are changed
% tic
if(A(u0,v0)==1)
Mold = A2old;
M = Mold;
% update u0 column of M
M(u0,:) = A2old(u0,:) + A(v0,:).*(1-A(u0,:));
M(u0,u0) = 0;
% update v0 column of M
M(:,v0) = A2old(:,v0) + A(:,u0).*(1-A(:,v0));
M(v0,v0)=0;
M(u0,v0) = 0;
else
Mold = A2old;
M = Mold;
% update u0 column of M
M(u0,:) = A2old(u0,:) - A(v0,:).*(1-A(u0,:));
M(u0,u0) = 0;
% update v0 column of M
M(:,v0) = A2old(:,v0) - A(:,u0).*(1-A(:,v0));
M(v0,v0)=0;
M(u0,v0) = A(u0,:)*A(:,v0);
end
% check that the answer is correct
fprintf('Diff between M and A2new = %f\n',norm(M - A2new,'fro'));


%% STEP 6 
% update the payoff, only updating those parts that changed
% Note when dividing by M(u,v), I do this:
% (M(u,v)~=0)/(M(u,v) + (M(u,v)==0)) 
% This is just a quick way to avoid a division by zero when M(u,v)==0
payoff = pold;
tic
if(A(u0,v0)==1)
hnew = (1 - A(u0,:)).*(M(u0,:)~=0)./(M(u0,:) + (M(u0,:)==0));
hold = (1 - A(u0,:)).*(Mold(u0,:)~=0)./(Mold(u0,:) + (Mold(u0,:)==0));
payoff = payoff + A(u0,:)'.*A*(hnew-hold)';
    
% payoff=pold;
% gnew = zeros(N,1);
% gold = zeros(N,1);
% for w=1:N
%     for v=1:N
%         gnew(w) = gnew(w) + A(u0,w)*A(w,v)*(1-A(u0,v))*...
%             (M(u0,v)~=0)/(M(u0,v) + (M(u0,v)==0));
%         gold(w) = gold(w) + A(u0,w)*A(w,v)*(1-A(u0,v))*...
%             (Mold(u0,v)~=0)/(Mold(u0,v) + (Mold(u0,v)==0));
%     end
%     payoff(w) = payoff(w) + (gnew(w) - gold(w));
% end
%norm(payofftest-payoff)
hnew = (1 - A(:,v0)).*(M(:,v0)~=0)./(M(:,v0) + (M(:,v0)==0));
hold = (1 - A(:,v0)).*(Mold(:,v0)~=0)./(Mold(:,v0) + (Mold(:,v0)==0));

payoff = payoff + A(:,v0).*A'*(hnew-hold);
payoff = payoff - A(u0,:)'.*A(:,v0)*...
       (Mold(u0,v0)~=0)/(Mold(u0,v0) + (Mold(u0,v0)==0)); 

% gnew = zeros(N,1);
% gold = zeros(N,1);
% for w=1:N
%     for u=1:N
%         gnew(w) = gnew(w) + A(u,w)*A(w,v0)*(1 - A(u,v0))*...
%             (M(u,v0)~=0)/(M(u,v0) + (M(u,v0)==0));
%         gold(w) = gold(w) + A(u,w)*A(w,v0)*(1 - A(u,v0))*...
%             (Mold(u,v0)~=0)/(Mold(u,v0) + (Mold(u,v0)==0));
%     end
%     payoff(w) = payoff(w) + (gnew(w) - gold(w));
%     payoff(w) = payoff(w) - A(u0,w)*A(w,v0)*...
%        (Mold(u0,v0)~=0)/(Mold(u0,v0) + (Mold(u0,v0)==0));    
% end
% norm(payofftest-payoff)

payoff(u0) = pold(u0);
payoff(u0) = payoff(u0)+ sum(A(:,u0).*(1 - A(:,v0)).*...
        (M(:,v0)~=0)./(M(:,v0) + (M(:,v0)==0)),1);
% for u=1:N
%     payoff(u0) = payoff(u0)+ A(u,u0)*(1 - A(u,v0))*...
%         (M(u,v0)~=0)/(M(u,v0) + (M(u,v0)==0));
% end
payoff(v0) = pold(v0);
payoff(v0) = payoff(v0)+ sum(A(v0,:).*(1 - A(u0,:)).*...
        (M(u0,:)~=0)./(M(u0,:) + (M(u0,:)==0)),2);
% for v=1:N
%     payoff(v0) = payoff(v0)+ A(v0,v)*(1 - A(u0,v))*...
%         (M(u0,v)~=0)/(M(u0,v) + (M(u0,v)==0));
% end


else
hnew = (1 - A(u0,:)).*(M(u0,:)~=0)./(M(u0,:) + (M(u0,:)==0));
hold = (1 - A(u0,:)).*(Mold(u0,:)~=0)./(Mold(u0,:) + (Mold(u0,:)==0));
payoff = payoff + A(u0,:)'.*A*(hnew-hold)';
    
% payoff=pold;
% gnew = zeros(N,1);
% gold = zeros(N,1);
% for w=1:N
%     for v=1:N
%         gnew(w) = gnew(w) + A(u0,w)*A(w,v)*(1-A(u0,v))*...
%             (M(u0,v)~=0)/(M(u0,v) + (M(u0,v)==0));
%         gold(w) = gold(w) + A(u0,w)*A(w,v)*(1-A(u0,v))*...
%             (Mold(u0,v)~=0)/(Mold(u0,v) + (Mold(u0,v)==0));
%     end
%     payoff(w) = payoff(w) + (gnew(w) - gold(w));
% end
%norm(payofftest-payoff)
hnew = (1 - A(:,v0)).*(M(:,v0)~=0)./(M(:,v0) + (M(:,v0)==0));
hold = (1 - A(:,v0)).*(Mold(:,v0)~=0)./(Mold(:,v0) + (Mold(:,v0)==0));

payoff = payoff + A(:,v0).*A'*(hnew-hold);
payoff = payoff + A(u0,:)'.*A(:,v0)*...
       (Mold(u0,v0)~=0)/(Mold(u0,v0) + (Mold(u0,v0)==0)); 

% gnew = zeros(N,1);
% gold = zeros(N,1);
% for w=1:N
%     for u=1:N
%         gnew(w) = gnew(w) + A(u,w)*A(w,v0)*(1 - A(u,v0))*...
%             (M(u,v0)~=0)/(M(u,v0) + (M(u,v0)==0));
%         gold(w) = gold(w) + A(u,w)*A(w,v0)*(1 - A(u,v0))*...
%             (Mold(u,v0)~=0)/(Mold(u,v0) + (Mold(u,v0)==0));
%     end
%     payoff(w) = payoff(w) + (gnew(w) - gold(w));
%     payoff(w) = payoff(w) - A(u0,w)*A(w,v0)*...
%        (Mold(u0,v0)~=0)/(Mold(u0,v0) + (Mold(u0,v0)==0));    
% end
% norm(payofftest-payoff)

payoff(u0) = pold(u0);
payoff(u0) = payoff(u0)- sum(A(:,u0).*(1 - A(:,v0)).*...
        (M(:,v0)~=0)./(M(:,v0) + (M(:,v0)==0)),1);
% for u=1:N
%     payoff(u0) = payoff(u0)+ A(u,u0)*(1 - A(u,v0))*...
%         (M(u,v0)~=0)/(M(u,v0) + (M(u,v0)==0));
% end
payoff(v0) = pold(v0);
payoff(v0) = payoff(v0)- sum(A(v0,:).*(1 - A(u0,:)).*...
        (M(u0,:)~=0)./(M(u0,:) + (M(u0,:)==0)),2);
% for v=1:N
%     payoff(v0) = payoff(v0)+ A(v0,v)*(1 - A(u0,v))*...
%         (M(u0,v)~=0)/(M(u0,v) + (M(u0,v)==0));
% end    
 
end
time_efficient(t) = toc;
edges(t) = sum(sum(A));
fprintf('time taken %d\n', time_efficient(t));

a = norm(payoff - pnew);
fprintf('Diff between payoff and pnew = %f\n',a);

payoff = payoff * 2 * alpha0;
% direct payoff
if A(u0,v0) == 1
    payoff(u0) = payoff(u0) + alpha0 - cost;
    payoff(v0) = payoff(v0) + alpha0 - cost;
else
    payoff(u0) = payoff(u0) - alpha0 + cost;
    payoff(v0) = payoff(v0) - alpha0 + cost;
 end
     % if the payoff decreases, change it back
if payoff(u0) < 0 || payoff(v0) < 0
    M = Mold;
    payoff = pold;
    if A(u0,v0) == 1
        A(u0,v0) = 0;
    else
        A(u0,v0) = 1;
    end
%     else
%         if A(u0,v0) == 1
%             A(v0,u0) = 1;          
%         else 
%             A(v0,u0) = 0;
%          end
end

end

function [payoff,time_efficient] = efficientPayoff2(A)
% test the payoff calculation by calculating the payoff of
% a graph A by adding one edge at a time and updating the 
% payoff after each addition of an edge
N = size(A,1);

% gather all the edges
[i,j]=find(A);
numEdges=length(i);
fprintf('Total number of edges in the graph is %d\n', numEdges);
maxEdges=input('Enter number of edges to test\n');



Mold = sparse([],[],[],N,N);
pold = zeros(N,1);

% empty A
A = sparse([],[],[],N,N);

%% STEP 3 
% Add one edge at a time and update the payoff for each addition
tic;
for e=1:min(numEdges,maxEdges)
    u0=i(e);
    v0=j(e);
    if (mod(e,1)==0)
        fprintf('%d\t%e\t%f\n', e,sum(pold),toc);
    end
        
    A(u0,v0)=1;
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    % Now update A2, only updating the rows and columns that are changed
    M = Mold;
    % update u0 column of M
    M(u0,:) = Mold(u0,:) + A(v0,:).*(1-A(u0,:));
    M(u0,u0) = 0;
    % update v0 column of M
    M(:,v0) = Mold(:,v0) + A(:,u0).*(1-A(:,v0));
    M(v0,v0) = 0;
    M(u0,v0) = 0;
    
    
    %%
    % update the payoff, only updating those parts that changed
    % Note when dividing by M(u,v), I do this:
    % (M(u,v)~=0)/(M(u,v) + (M(u,v)==0))
    % This is just a quick way to avoid a division by zero when M(u,v)==0
    payoff = pold;

    hnew = (1 - A(u0,:)).*(M(u0,:)~=0)./(M(u0,:) + (M(u0,:)==0));
    hold = (1 - A(u0,:)).*(Mold(u0,:)~=0)./(Mold(u0,:) + (Mold(u0,:)==0));
    payoff = payoff + A(u0,:)'.*A*(hnew-hold)';
  
    hnew = (1 - A(:,v0)).*(M(:,v0)~=0)./(M(:,v0) + (M(:,v0)==0));
    hold = (1 - A(:,v0)).*(Mold(:,v0)~=0)./(Mold(:,v0) + (Mold(:,v0)==0));
    payoff = payoff + A(:,v0).*A'*(hnew-hold);
    
    payoff = payoff - A(u0,:)'.*A(:,v0)*...
        (Mold(u0,v0)~=0)/(Mold(u0,v0) + (Mold(u0,v0)==0));
    
    payoff(u0) = pold(u0);
    payoff(u0) = payoff(u0)+ sum(A(:,u0).*(1 - A(:,v0)).*...
        (M(:,v0)~=0)./(M(:,v0) + (M(:,v0)==0)),1);
    
    payoff(v0) = pold(v0);
    payoff(v0) = payoff(v0)+ sum(A(v0,:).*(1 - A(u0,:)).*...
        (M(u0,:)~=0)./(M(u0,:) + (M(u0,:)==0)),2);
    
    
    pold = payoff;
    Mold = M;
end
%% STEP 4 
% Finally check if the payoff calculated edge-by-edge
% equals the payoff calculated on the full graph
A2=((A*A).*(A==0));
A2 = A2 - sparse(1:N,1:N,diag(A2));
% Secondly calculate the payoff
G = (1-A).*(A2~=0)./(A2+(A2==0));
H = A*G';
fullpayoff = sum(A'.*H,2);

fprintf('Difference between delta payoff and full payoff=%f\n', ...
    norm(fullpayoff-payoff));

time_efficient = toc;
end

function [payoff,time_inefficient] = inefficientPayoff2(A)
% test the payoff calculation by calculating the payoff of
% a graph A by adding one edge at a time and updating the 
% payoff after each addition of an edge
N = size(A,1);

[i,j]=find(A);
numEdges=length(i);
fprintf('Total number of edges in the graph is %d\n', numEdges);

Mold = sparse([],[],[],N,N);
pold = zeros(N,1);

A = sparse([],[],[],N,N);

% Add one edge at a time and update the payoff for each addition
tic;
for e=1:1000
    u0=i(e);
    v0=j(e);
    if (mod(e,1)==0)
        fprintf('%d\t%e\t%f\n', e,sum(pold),toc);
    end
    if (A(u0,v0)==1)
        fprintf('adding to existing edge\n');
        exit;
    end
    A(u0,v0)=1;
    % Firstly calculate A2
    A2=((A*A).*(A==0));
    A2 = A2 - sparse(1:N,1:N,diag(A2));
    
    % Secondly calculate the payoff
    
    G = (1-A).*(A2~=0)./(A2+(A2==0));
    H = A*G';
    payoff = sum(A'.*H,2);
    pold = full(payoff);
end
time_inefficient = toc;
end
