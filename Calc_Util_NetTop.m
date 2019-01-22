[Adj_, nd_coord_, Net_Top, forward_weights, backward_weights, delta_t]= Network_Topology_Generation(); 
nodes = linspace(1,10,10); 
node_subsets = PowerSet(nodes); 

reliabilities = zeros(1,length(node_subsets)); 
for u = 2: length(node_subsets) %start at 2 to avoid empty set 
    %disp(node_subsets{1,u}); 
    val = throughput(Net_Top, node_subsets{1,u}, delta_t);
    fprintf("COALITION: "); 
    disp(node_subsets{1,u}); 
    fprintf("THROUGHPUT: %d\n\n", val); 
    %disp(size(reliabilities)); 
    reliabilities(u) = val; 
end

% How do we fairly distribute the gains among all the nodes. Whay payoff
% can nodes reasonably expect from cooperation? Shapley value is one way to
% distribute the total gains to players. 

%x_i(v) = sum over all subsets not including i. 
% utility distribution within coalitions 

node_payoff = zeros(1, length(nodes)); 
n = length(nodes); 
n_fac = factorial(n);

for h = 1: 10
    x_v = 0; 
    for g = 1: length(node_subsets)
        curr_subset = node_subsets{1,g};  
        found = find(curr_subset == h);
  
        if isempty(found) 
            S_mag = length(curr_subset); 
            utility_S = throughput(Net_Top, curr_subset, delta_t); 
            utility_U = throughput(Net_Top, union(curr_subset, h), delta_t); 
        
            x_v = x_v + (factorial(S_mag)*factorial(n - S_mag - 1)*(utility_U - utility_S))/n_fac;
        end
    end
    node_payoff(1,h) = x_v; 
end

disp(node_payoff); 
    

%generating powerset of N 
function [ P ] = PowerSet(S)
    n = numel(S); 
    x = 1:n; 
    P = cell(1,2^n); 
    p_ix = 2; 
    for nn = 1:n 
        a = combnk(x,nn); 
        for j=1:size(a,1) 
            P{p_ix} = S(a(j,:)); 
            p_ix = p_ix + 1; 
        end
    end
end

function [reliability] = throughput(N, S, delta_t) 
% N = grand coalition --> entire graph 
% S is a subset of the nodes in N
% Q_ab is the number of packets required to be sent from a source node
% to a destination node 

coalition = subgraph(N, S);
C_adj = adjacency(coalition); 

% for each source, destination pair in coalition
reliability = 0; 
for i = 1: length(S)
    for j = 1: length(S)
        if i ~= j
            %generating new Q_ab for each time an edge is calculated is
            %wrong -- we should have a unique, stable Q_ab for each edge 
            if findedge(N,S(i),S(j)) ~= 0
                if i < j
                    Q_ab = N.Edges.Q_abf(findedge(N,S(i),S(j))); 
                elseif j < i 
                    Q_ab = N.Edges.Q_abb(findedge(N,S(i),S(j)));
                end
            else
                Q_ab = 0; 
            end 
            
            [max_k] = t_k(coalition, C_adj, i, j); 
            reliability = reliability + (Q_ab*max_k);
        end
    end
end

reliability = reliability/delta_t; 
end 

function [most_reliable_path] = t_k(sub_graph, Adj_, src, dest) 
%takes a subgraph of Net_Top and outputs the score of the most reliable path: max_path
%We have a subset S of a grand coalition N. Given the source destination
%pair in Subgraph, we want to find the "cost" (reliability) of the most
%secure path

%1) get all paths between src and dest in the subgraph
%2) calculate reliability of each path 
%3) get the max of the reliability
%4) Note that the index of the max reliability corresponds to the index of
%the path with the max reliability 

[paths] = pathbetweennodes(Adj_, src, dest);  
path_scores = zeros(1,length(paths)); 

if length(paths) >= 1 
    % for each path, multiply all the edgeweights 
    % we want to return the largest multiplication of edgeweights 
    for i = 1: length(paths) %loop through all paths
        path_score = 1;
        path = paths{i,1};
        for j = 1: length(path) - 1
            if path(j) < path(j+1)
                % path contains the nodes along path from src to dest. We go to
                % each edge in the path and get either its forward or backward
                % weight from the subgraph
                path_score = path_score*sub_graph.Edges.Indf(findedge(sub_graph,path(j),path(j+1))); 
            else
                path_score = path_score*sub_graph.Edges.Indb(findedge(sub_graph,path(j),path(j+1)));
        
            end
            
        end
        path_scores(1,i) = path_score;
    end
    %disp(path_scores);
    most_reliable_path = max(path_scores); 
else
    most_reliable_path = 0;  
end 
end 



