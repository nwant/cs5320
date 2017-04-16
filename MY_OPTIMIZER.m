% jumping swarm algorithm
%----------------------------------------------
function MY_OPTIMIZER(FUN, DIM, ftarget, maxfunevals)
	% Set algorithm parameters
	popsize = 40;
	c1 = 1.4944;%2;
	c2 = 1.4944;%2;
	w = 0.792;
	xbound = 5;
	%vbound = 5;
	% Allocate memory and initialize
	xmin = -xbound * ones(1,DIM);
	xmax = xbound * ones(1,DIM);
	%vmin = -vbound * ones(1,DIM);
	%vmax = vbound * ones(1,DIM);
	%A = []																% no members of pop better than gbest yet	
	x = 2 * xbound * rand(popsize,DIM) - xbound;  % initial population for class B
	%v = 2 * vbound * rand(popsize,DIM) - vbound;
	%pbest = x;
	% update pbest and gbest
	% 
	cost_x = feval(FUN, x');
	%disp(cost_p);
	%disp("");
	[cost,index] = min(cost_x);
	cost_g = cost;																	 % keep track of the 
	gbest = x(index,:);															 % assign member of pop with best fitness as gbest
	%disp(gbest);
	garc = [gbest] 																	 % gbest archive
	disp(garc);
	maxfunevals = min(1e5 * DIM, maxfunevals);
	maxiterations = ceil(maxfunevals/popsize);
	for iter = 2 : maxiterations	

		% jump toward gbest as a swarm
		x = x + rand(popsize,DIM).* (repmat(gbest,popsize,1)-x);

		% Clamp position - Absorbing boundaries
		% Set x to the boundary
		s = x < repmat(xmin,popsize,1);
		x = (1-s).*x + s.*repmat(xmin,popsize,1);
		b = x > repmat(xmax,popsize,1);
		x = (1-b).*x + b.*repmat(xmax,popsize,1);

		% Update pbest and gbest if necessary
		cost_x = feval(FUN, x');
		%s = cost_x<cost_p;
		%cost_p = (1-s).*cost_p + s.*cost_x;
		%s = repmat(s',1,DIM);
		%pbest = (1-s).*pbest + s.*x;
		[cost,index] = min(cost_x);
		
		% has one or more solutions found a new gbest?
		if cost < gbest
			new_gbest = x(index,:);
			% select a random gbest from the archive
			% yes; prepare to relocate all solutions that have matched or beat gbest this iteration
			idx = randi(length(garc), popsize, 1); 				% generate a potential gbest archive index for each member of the population
			for i = 1 : popsize
				if cost_x(i) >= gbest:  										% for each member of the pop that has matched or beaten gbest this iteration:
					x(i) = mutate(gbest, gargc(idx(i)))				% create a new mutated "gbest" using this gbest and the archive and reassign position to this member.
				end
			end				
				
			% update gbest to new gbest and archive old 
		 	garc = [garc; gbest];
			gbest = new_gbest; 
		end

		% Exit if target is reached
		if feval(FUN, 'fbest') < ftarget
			break;
		end
	end
end

function mutate(g1, g2)
	rand_indices = randsample(1:length(G), length(x));
	grand = G(rand_indices,:);
end
