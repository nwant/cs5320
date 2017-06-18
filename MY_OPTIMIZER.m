% Jumping Swarm Algorithm  (JSO)
%=============================================
%	by Sulejman Uvalic and Nathan Want 
% CS5320-G01
% Spring 2017
%----------------------------------------------
%**
% MY_OPTIMIZER (JSO)
%-----
% The JSO algorithm created to work with COCO framework.
%
% Inputs:
% FUN..........[function] the fitness function (provided by COCO)
% DIM..........[number] the number of dimensions 
%              (provided by COCO, defined in exampleexperiment.m)
% ftarget......[number] the target fitness (provided by COCO)
% maxfunevalus.[number] the maximum number of function evaluations alotted 
%							 (provided by COCO, defined in exampleexperiment.m)
%
function MY_OPTIMIZER(FUN, DIM, ftarget, maxfunevals)
	% Set algorithm parameters
	swarmsize = 4000;	% total swarm size
	xbound = 5;			% [-xbound, +xbound] for each dimension, as defined by coco
	garc = []; 			% gbest archive (all positions that were gbest at one point in time) 

	% generate the initial swarm
	% initial swarm will be a [swarmsize x DIM] matrix
	% each member of the swarm will be 
	xmin = -xbound * ones(1,DIM);
	xmax = xbound * ones(1,DIM);
	x = 2 * xbound * rand(swarmsize,DIM) - xbound;  % set all dim for each member [-xbound, xbound]

	% Determine gbest and start gbest archive
	cost_x = feval(FUN, x');									 % evaluate the fitness of each swarm member
	[bestcost,index] = min(cost_x);						 % determine gbest fitness (bestcost)
	gbest = x(index,:);												 % assign member of swarm with best fitness as gbest

	% set termination criteria
	maxfunevals = min(1e8 * DIM, maxfunevals); % the maximum evaluations of fitness alotted
	maxiterations = ceil(maxfunevals/swarmsize); % the maximum "jumps" allotted	

	% perform algorithm until the optimal solution is found 
	% or we max out the number of allotted iterations. 
	for iter = 2 : maxiterations							 	

		% jump to new position and "clamp" any positions that fell outside out bounds.
		x = jump(x, swarmsize, gbest, DIM); 
		x = clamp(x, swarmsize, xmin, xmax);

		% determine the fitness of new position for each member
		% and the best fitness from within the swarm's new jump
		cost_x = feval(FUN, x');
		[cost,index] = min(cost_x);

		% is at least one position more fit than the current gbest?	
		if cost < bestcost  
			% yes; archive the previous gbest
			new_gbest = x(index,:); 	
		 	garc = [garc; gbest];			 
																										
			% Relocate all solutions that have matched or beat gbest this iteration
			x = relocate_best(x, swarmsize, garc, bestcost, new_gbest, cost_x);

			% set the new gbest for next jump
			gbest = new_gbest; 
			bestcost = cost;
		end

		% have we reached the target (criteria determined by COCO framework)
		if feval(FUN, 'fbest') < ftarget
			break;
		end
	end
end


%**
% jump 
%-----
% generate a new swarm by performing a "jump" toward the provided gbest.
%
% Inputs:
%	x..........[matrix; swarmsize x DIM] the current swarm
% swarmsize..[number] the size of the current swarm
% gbest......[matrix; 1 x DIM] the current gbest
% DIM........[number] the number of dimensions for each solution (member in swarm)
%
% Returns: [matrix; swarmsize x DIM] the swarm after the jump
%
function x = jump(x, swarmsize, gbest, DIM)
	% for each swarm member, somewhere between its current position and the gbest position.	
	x = x + rand(swarmsize,DIM).* (repmat(gbest,swarmsize,1)-x);
end


%**
% clamp  
%-----
% clamp a swarm (e.g. absorbing the boundaries) as provided when calling this function.
%
% Inputs:
%	x..........[matrix; swarmsize x DIM] the swarm to be clamped
% swarmsize..[number] the size of provided swarm
% xmin.......[matrix; 1 x DIM] matrix with all cells equal to the minimum dim size
% xmax.......[matrix; 1 x DIM] matrix with all cells equal to the maximum dim size
%
% Returns: [matrix; swarmsize x DIM] the swarm after clamping
%
function x = clamp(x, swarmsize, xmin, xmax)
	% clamp any dimensions for solutions that have dimensions more negative than the min dim size
	s = x < repmat(xmin,swarmsize,1);						
	x = (1-s).*x + s.*repmat(xmin,swarmsize,1); % 
	% clamp any dimensions for solutions that have dimensions more positive than the max dim size
	b = x > repmat(xmax,swarmsize,1);
	x = (1-b).*x + b.*repmat(xmax,swarmsize,1);
end


%**
% relocate_best  
%-----
% classify each swarm member by determining if it has matched or beat the provided gbest
% in terms of fitness, or it has not beat or matched the provided gbest. For all members
% that have matched or beat the current gbest, reposition them by assigning a position
% derived from the mutation of the current gbest and one gbests in the provided gbest
% archive.
%
% Inputs:
%	x..........[matrix; swarmsize x DIM] the swarm to be clamped
% swarmsize..[number] the size of provided swarm
% xmin.......[matrix; swarmsize x DIM] matrix with all cells equal to the minimum dim size
% xmax.......[matrix; swarmsize x DIM] matrix with all cells equal to the maximum dim size
% garc.......[matrix; # of archived gbests x DIM] the gbest archive
% bestcost...[number] the fitness of the current gbest
% gbest......[matrix; 1 x DIM] the position that has the best fitness found thus far.
% cost_x.....[matrix; swarmsize x 1] the fitness for each member in the swarm
%
% Returns: [matrix; swarmsize x DIM] the swarm after relocating the best swarm members
%
function x = relocate_best(x, swarmsize, garc, bestcost, gbest, cost_x)
	[rows, col] = size(garc);
	% generate a potential gbest archive index for mutation for each member of the swarm				
	idx = randi(rows, swarmsize, 1); 
	% for each member of the swarm that has matched or beaten gbest this iteration:
	for i = 1 : swarmsize
		% is this member part of the "best" class?
		if cost_x(i) <= bestcost  										  
			% yes; relocate to a mutated position based on current gbest and a random gbest from the archive.
			m = mutate(gbest, garc(idx(i), :));				
			x(i,:) = m;
		end
	end				
end


%**
% mutate
%-----
% Create a mutant position by using attributes of the provided positions.
% Attributes in the first provided postion are favored over the attributes 
% of the second provided solution.
%
% Inputs:
% x1....[matrix; 1 x DIM] the position to favor the attributes for when 
%       creating the mutant position
% x2....[matrix; 1 x DIM] another position whose attributes will also be 
%       utilized when creating the mutant position.
%
% Returns: [matrix; 1 x DIM] the mutant position
function m = mutate(x1, x2)
	% randomly select an attribute from x2
	idx = randi(length(x1));
	% swap this attribute with the same attribute in x1
	x1(idx) = x2(idx);
	% return the mutated x1.
	m = x1;
end
