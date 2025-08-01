This model introduces loose priors on selex parameters in an attempt to stabilize model estimates.

	--Priors are normal with mean of ln(1) and sd of 5.

Priors on trawl fishery selex are slightly tighter, since the gamma seems a bit more unstable and only length data are available (became very unstable once age comps split).

	--gamma: mean of ln(2), sd of 1.
	--bmax: mean ln(1) and sd of 2.0

The original linking of selex pars from the Cont model are also removed and replaced with three linkages:
	--Fixed gear fishery early period (IFQ pre-1995) selex deltas are linked between males and females (1 delta estimated). 
	--JPN LLS selex deltas are linked between males and females (1 delta estimated).
	-- Trawl fishery gamma function bmax pars are linked between males and females (1 bmax estimated).
These linkages were chosen a priori because there are extremely limited comp data from which to estimate selex pars from the JPN LLS and from the fishery before 1995, and trawl fishery
selex is highly correlated with M.