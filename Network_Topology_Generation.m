function [A, nd_coord, W_T, Indf_, Indb_, delta_t]= Network_Topology_Generation(radio_range,xmax,xmin,ymax,ymin,npoints, alpha,delta_t)

% Inputs: 
%   domain - bounds for the region.
%   radio_range - max distance between nodes for them to be able to comm
%   npoints - number of points to be randomly distributed 
%
% Outputs:
%   adj_matr - adjacency matrix of the graph of the topology
%   nd_coord - coordinates of the nodes

%Generating Initial Wireless Topology 
%Nodes within range of one another are connected by an edge 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if (nargin<4) % default parameter values
   radio_range = 32; 
   xmax = 100;
   xmin = 0;
   ymax = 100;
   ymin = 0; % bounds for the "geographical" domain
   npoints = 10; 
   alpha = .6; %weight of direct experience of communication between nodes. Weight of indirect experience is given by 1-alpha
   delta_t = 1000; 
 end

 % given the number of points, nodes are uniformly distributed
 nd_coord = rand(npoints, 2); % 2 cols for (x,y) coordinates
 nd_coord(:, 1) = nd_coord(:, 1)*(xmax-xmin)+xmin;
 nd_coord(:, 2) = nd_coord(:, 2)*(ymax-ymin)+ymin;
 
 % now that we have our points randomly distributed, we want to create
 % an adjacency matrix that connects points if they are in "range" with 
 % one another 
 
 % we can plot the network topology on a euclidean plane using an adjacency
 % matrix, and represent the topology using a graph. For simplicity's sake,
 % we assume that the graph is undirected; that is, that if node_1 can
 % communicated with node_2, then node_2 can communicated with node_1. We
 % have two Trust probability weights because each source-destination pair
 % (s,d) has a trust value, and T(s_1,d_1) != T(d_1,s_1). 
 
 % create a matrix with all possible distances
 x_rep = repmat(nd_coord(:, 1), 1, npoints); 
 y_rep = repmat(nd_coord(:, 2), 1, npoints); 
 dist_matr = sparse(triu(((x_rep-x_rep').^2 +(y_rep-y_rep').^2).^0.5, 1)); 
 
 %generate the adjacency matrix 
 A = (dist_matr <= radio_range) & (dist_matr > 0); 
 A = A + A'; 
 
 % creating undirected graph 
 W_T = graph(A); % we assume that each edge is bidirectional 
 trust1_ = .5*rand(W_T.numedges,1)+.2; % 
 trust2_ = .5*rand(W_T.numedges,1)+.2; % Some problems here, because ideally this reliability measure would be a function of distance  
 W_T.Edges.Trust1 = trust1_; %"forward edges" direct experience
 W_T.Edges.Trust2 = trust2_; %"backward edges" direct experience 
 % disp(W_T.Edges); 
 % plot 

 source_ = linspace(1,10,10);
 Dist_ = zeros(W_T.numedges,1); 
%  zero_edges = zeros(20,6);
%  zero_edge_count = 1; 
 edge_number = 1; 
 for i = 1: length(source_) 
     %Succ = successors(W_T,i);
     N = neighbors(W_T,i); 
     if ~isempty(N) %~isempty(Succ)
         %indirect_exp = 0; 
         
         %for each neighbor of i
         for j = 1: length(N)
             %don't check any nodes < i in N because those edges have
             %already been checked 
             %dest = N(j); 
             if(N(j) > i)
                 
                 %it appears that occassionally Dist = 0. Let's record for
                 %which edges that seems to be the case and look at the
                 %coordinates of the nodes joined by the edge
%                  if (((nd_coord(i,1)-nd_coord(j,1))^2 +(nd_coord(i,2)-nd_coord(j,2))^2)^0.5) == 0
%                      zero_edges(zero_edge_count,1) = i; 
%                      zero_edges(zero_edge_count,2) = j; 
%                      zero_edges(zero_edge_count,3) = nd_coord(i,1);
%                      zero_edges(zero_edge_count,4) = nd_coord(i,2); 
%                      zero_edges(zero_edge_count,5) = nd_coord(j,1); 
%                      zero_edges(zero_edge_count,6) = nd_coord(j,2); 
%                  end
                 Dist_(edge_number,1) = (((nd_coord(i,1)-nd_coord(N(j),1))^2 +(nd_coord(i,2)-nd_coord(N(j),2))^2)^0.5);
                 edge_number = edge_number + 1; 
             end 
             
         end
     end 
 end
 W_T.Edges.Dist = Dist_; 
 
 disp(W_T.Edges); 
 
 Indf_ = zeros(W_T.numedges,3); 
 Indb_ = zeros(W_T.numedges,3); %  
 forward_count = 1; 
 backward_count = 1;
 for i = 1: length(source_)
     
     N = neighbors(W_T,i); 
     
     if ~isempty(N)
         
         for j = 1: length(N)
             %N(j) becomes the destination node, so the edge 
             %we are constructing P_ij for (i,N(j))
             %N(j) = j
             if(N(j) > i) %forward edges, N(j) = j (destination)
                 P_ij_forward = 0; 
                 num_fedges = 0;
                 for k = 1 : length(N) % N(k) = l (neighbors of source)
                     if N(k) ~= N(j) 
                         P_il_f = W_T.Edges.Trust1(findedge(W_T,i,N(k)));
                         P_li_f = W_T.Edges.Trust2(findedge(W_T,i,N(k))); 
                         P_lj_f = 1; 
                         if findedge(W_T,N(k),N(j)) ~= 0 %p_lj will be undefined if there is no edge between N(k) and N(j), so we have to initialize P_lj_f = 1 
                            if N(k) < N(j) 
%                                 fprintf("N(K)_f = %d\nN(J)_f = %d\n", N(k), N(j));
%                                 fprintf("edge index = %d\n",findedge(W_T,N(k),N(j)));
                                P_lj_f = W_T.Edges.Trust1(findedge(W_T,N(k),N(j)));
                            else
%                                 fprintf("N(K)_b = %d\nN(J)_b = %d\n", N(k), N(j)); 
%                                 fprintf("edge index = %d\n",findedge(W_T,N(k),N(j)));
                                P_lj_f = W_T.Edges.Trust2(findedge(W_T,N(k),N(j)));
                            end
                            P_ij_forward = P_ij_forward + (P_il_f*P_li_f*P_lj_f); 
                            num_fedges = num_fedges + 1;
                         end
                     end
                 end
                 % if P_ij_forward remains 0, then there is no indirect
                 % experience to consider
                 P_ij_forward = P_ij_forward/length(N); %dividing by zero sometimes
                 Indf_(forward_count,1) = i; %source
                 Indf_(forward_count,2) = N(j); %dest
                 Indf_(forward_count,3) = P_ij_forward; 
                 
                 forward_count = forward_count + 1; 
                 
             else %backwards edges, again N(j) = dest, N(l) = neighbor of i 
                 P_ij_backward = 0; 
                 num_bedges = 0; 
                 for l = 1 : length(N) % N(k) = l (neighbors of source)
                     
                     if N(l) ~= N(j) 
                         P_il_b = W_T.Edges.Trust2(findedge(W_T,i,N(l)));
                         P_li_b = W_T.Edges.Trust1(findedge(W_T,i,N(l))); 
                         P_lj_b = 1; 
                         %we are looking at the edge (N(l),N(j))
                         %if the edge exists, we check whether it is a
                         %forward or backwards edge. Forward if N(l) < N(j)
                         %and backwards if N(j) < N(l) 
                         if findedge(W_T,N(l),N(j)) ~= 0
                            if N(l) < N(j) 
%                                 fprintf("N(K)_f = %d\nN(J)_f = %d\n", N(l), N(j));
%                                 fprintf("edge index = %d\n",findedge(W_T,N(l),N(j)));
                                P_lj_b = W_T.Edges.Trust2(findedge(W_T,N(l),N(j)));
                            else
%                                 fprintf("N(K)_b = %d\nN(J)_b = %d\n", N(l), N(j)); 
%                                 fprintf("edge index = %d\n",findedge(W_T,N(l),N(j)));
                                P_lj_b = W_T.Edges.Trust1(findedge(W_T,N(l),N(j)));
                            end
                         end
                         
                         P_ij_backward = P_ij_backward + (P_il_b*P_li_b*P_lj_b); 
                         num_bedges = num_bedges + 1;
                     end
                 end
                 
                 P_ij_backward = P_ij_backward/length(N); 
                 Indb_(backward_count,1) = i; %source
                 Indb_(backward_count,2) = N(j); %dest
                 Indb_(backward_count,3) = P_ij_backward; 
                 backward_count = backward_count+1; 
             end
         end
     end
 end

 %sort Indb_  
 Indb_ = sortrows(Indb_,2); 
 
%  disp(W_T.numedges); 
%  disp(size(Indf_));
%  disp(size(Indb_)); 
%  disp(Indf_);
%  disp(Indb_); 

%  disp(zero_edges); 
 
 W_T.Edges.Indf = Indf_(:,3); 
 W_T.Edges.Indb = Indb_(:,3);

 %We have P_ij for forward edges
 %We have D for all edges
 %We can now calculate the desired weight of edges which is P_ij/D^2. Note
 %however that this must be done for forward edges and for backward edges. 
 %We therefore will have both Weightf, and Weightb 

 Weightf_ = (alpha.*W_T.Edges.Trust1 + (1-alpha).*W_T.Edges.Indf)./W_T.Edges.Dist.^2; 
 Weightb_ = (alpha.*W_T.Edges.Trust2 + (1-alpha).*W_T.Edges.Indb)./W_T.Edges.Dist.^2;
 
 W_T.Edges.Weightf = Weightf_; 
 W_T.Edges.Weightb = Weightb_; 
 W_T.Edges.Weight = [];
 
 %associate each edge with a "required number of packets to be transmitted"
 %# of packets is taken to mean #packets/sec
 Q_ab_f = (600000)*rand(W_T.numedges,1)+200000;
 Q_ab_b = (600000)*rand(W_T.numedges,1)+200000;
 W_T.Edges.Q_abf = Q_ab_f; 
 W_T.Edges.Q_abb = Q_ab_b; 

 disp(W_T.Edges);
 
 figure(1);
 clf;
 hold on;
 plot(nd_coord(:, 1), nd_coord(:, 2), 'o',"linewidth",2,"MarkerSize", 2);
 gplot(A, nd_coord);
 title("Randomly Generated Network Topology", "fontsize", 15); 
 hold off;
 
 %Now we have P_ij, D_ij, delta_t, and a graph 
 %still confused about whether to use number of forward and back edges or
 %num_neighbors as |NB_i|. We are close to being able to define the
 %throughput function, and closer still to implementing simplex. 
 
 %Throughput function
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
 
 
 
 
 
