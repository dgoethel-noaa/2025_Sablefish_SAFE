Main changes:
--Use true likelihoods with variance weighting and removing lambdas for indices/priors/catch (still use for comps, eg Francis reweights), 
particularly for the catch (ADMB just used sum of squares with a lambda of 50)
--Set sigR to 0.9 because this more accurately reflects the ADMB SRR penalty function that used a lambda of 1.5 and sigR ~ 1.1 (inherently implying sigR was more like 0.9 than 1.1).
--Set sigC to 0.075 to reflect that catch is well known, but still some uncertainty in magnitude.
--Multiply all CPUE CV by 2 and survey indices CV by 2, which better reflects the true information content:
		-- CPUE CV was 0.1 now 0.2
		-- JPN LLS CV was 0.05 now 0.1
		-- LLS CV was scaled to ~0.05 now 0.1
		-- TS CV was scaled to ~0.135 now ~0.27
--Set Fprior SD to 1 (default)

Results allow model to fit catch on par with ADMB, but much better fits to the LLS index (with lower CVs, there was too much tension among data components leading to instability
and inability to actually fit the catch and LLS simulataneously).