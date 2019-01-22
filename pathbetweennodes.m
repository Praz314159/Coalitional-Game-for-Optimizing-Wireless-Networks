function [pth] = pathbetweennodes(adj, src, snk, verbose)


% major problem -- for some reason, all paths between nodes are not being
% provided. For example, fig 1: (1,3) only produces (1,3), even though
% (1,3) isn't a path. (1,6,2,3), (1,4,6,2,3), (1,7,4,6,2,3) these should
% all be paths. 

if nargin < 4
    verbose = false;
end

n = size(adj,1);

stack = src;

stop = false;

pth = cell(0);
cycles = cell(0);

next = cell(n,1);
for in = 1:n
    next{in} = find(adj(in,:));
end

visited = cell(0);

pred = src;
while 1
    
    visited = [visited; sprintf('%d,', stack)];
    
    [stack, pred] = addnode(stack, next, visited, pred);
    if verbose
        fprintf('%2d ', stack);
        fprintf('\n');
    end
    
    if isempty(stack)
        break;
    end
    
    if stack(end) == snk
        pth = [pth; {stack}];
        visited = [visited; sprintf('%d,', stack)];
        stack = popnode(stack);
    elseif length(unique(stack)) < length(stack)
        cycles = [cycles; {stack}];
        visited = [visited; sprintf('%d,', stack)];
        stack = popnode(stack);  
    end

end


function [stack, pred] = addnode(stack, next, visited, pred)

newnode = setdiff(next{stack(end)}, pred);
possible = arrayfun(@(x) sprintf('%d,', [stack x]), newnode, 'uni', 0);

isnew = ~ismember(possible, visited);

if any(isnew)
    idx = find(isnew, 1);
    stack = str2num(possible{idx});
    pred = stack(end-1);
else
    [stack, pred] = popnode(stack);
end


function [stack, pred] = popnode(stack)

stack = stack(1:end-1);
if length(stack) > 1
    pred = stack(end-1);
else
    pred = [];
end
