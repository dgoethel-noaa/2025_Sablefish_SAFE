This model fixes a bunch of legacy ADMB code that either was a due to coding bugs or was inconsistently applied across data sources.

	--The natural mortality male offset parameter (Mdev) was set to 0.0. This had been accidentally set to -0.008 in 2023 during a jitter analysis.
	--Lognormal additive constants in the catch likelihood were removed because catch was no longer fit if there was a 0 value (no possibility of log(0)). Unusually large values were 	being used (especially for the trawl survey ==0.8).
	--There were coding errors in how the fishery selectivity was being used to calculate the CPUE index in ADMB (only female selectivity not a combo of male and female). This was 	fixed.
	--The timing of the survey for some of the composition data was incorrectly set at 0 instead of mid-year. This was corrected.
	--The application of aging error and size-age transitions (for length data) were being applied before or after normalizing to sum to 1 across ages. A consistent approach is now 	utilized where aging or size-age transition is applied then normalization happens after (this is deemed good practice).