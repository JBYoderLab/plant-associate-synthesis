Population structure in plants and their associates
===================================================

Notes towards a review or synthesis


2019.06.24
----------

Startup --- blurted out this idea to Lynda Delph in a New Phytologist editorial meeting. Need to do the lit search to determine
- how many papers have pop gen structure (not phylogeography?) for a plant and at least one associate
- what kinds of patterns do we expect?
	- host/parasite theory
	- mutualist theory
	- &c
	
Starting points for the lit search
- citations of Smith et al 2010 *PLOS ONE*?
- citations of STRUCTURE that use multiple species (woof but)

If I can get enough multi-species papers, we could do some synthesis:
- correlation of Fst or the like?
- correlation of optimal k?


2019.06.26
----------

Lit-searching. (Just doing title/abstract checks; will review and cull later.)

Citing Smith *et al.* 2010, likely *data*
```
https://doi.org/10.1111/bij.12393
https://doi.org/10.1111/evo.12924
https://doi.org/10.1111/mec.12318
https://doi.org/10.1111/j.1365-294X.2012.05618.x
https://doi.org/10.1111/mec.12406
https://doi.org/10.1111/j.1365-294X.2012.05700.x
https://doi.org/10.1111/mec.13438
https://doi.org/10.1111/mec.14137
https://doi.org/10.1111/jbi.12456
https://doi.org/10.1111/jeb.12980
https://doi.org/10.3389/fevo.2019.00206
https://doi.org/10.1111/mec.15111 --- maybe? (fungal)
```

Citing Smith *et al.* 2010, *conceptual*
```
https://doi.org/10.1073/pnas.1604338113
https://doi.org/10.1093/cz/zow018
```


2019.06.28
----------

Reading through results from Wednesday, gonna use the following criteria
- population genetics analysis (STRUCTURE, Fst)
- of a plant and at least one associate

Good prototype is [Boyle *et al.* (2019, *Frontiers Ecol Evol*)](https://doi.org/10.3389/fevo.2019.00206).


2019.06.03
----------

Little more lit-searching: GSchol, terms `comparative population genetic host plant`

Prospective data. 

```
https://doi.org/10.1007/s10709-007-9153-6
https://doi.org/10.1111/j.0014-3820.2004.tb00457.x
https://doi.org/10.1046/j.0962-1083.2002.01462.x
https://doi.org/10.1111/j.1365-294X.2011.05296.x
```

For followup

```
https://doi.org/10.1098/rspb.2005.3363
https://doi.org/10.2307/2960598
```


2019.07.01
----------

Writing up proposal for LD; I have 8 papers so far with usable data, covering more species than that --- I am pretty sure I can expect that number to expand as I try new search terms sets. 


2019.10.30
----------

Have got about 16 papers with results from Alby, Cate and Mikhail's searching. Set up shared folder and data workbook with separate sheets for metadata, Fst distances, and K values. By way of starting data entry I'm going through the process of getting Fst for Joshua tree and its pollinators.


2020.02.07
----------

Have got what may be a final-ish dataset! Total of 15 papers, 15 plant taxa covered. 

- Still need to figure out geographic distances for the Evans (2013) dataset, but I'm prepared to take a whack at it with Google Maps, I think.
- Also, one cases in which the correlation between plant and associate Fst is significantly *negative*? WTF? Triponez. This is the one that came packaged with a script to calculate genetic distances ...

Thinking about next-level stuff: could I run BEDASSLE on this data? BEDASSLE requires standard Coop Lab input, allele-count data assuming biallelic loci, which I mostly don't have. However, BEDASSLE uses that data pretty early on to generate a pairwise matrix of Fst distances, and I *do* have those. If I can insert my own genetic distance matrices into the analysis, I could in principle make this work. Emailed Gideon Bradburd about it.

*later* GOT the distances from Evans 2013, thanks to some creative XYScanning. Have re-run all the data extraction from Triponez and ... the pattern persists. Ugh.


2020.02.08
----------

Went back to the text of Triponez 2014 and .. wow the negative plant-associate correlation is *reported in the paper*. Table 2 has "phylogeographical congruence" correlations negative for plant-M.europaea. I guess it's real! And with that, I can move forward with analysis.


2020.02.17
----------

First attempt to set up BEDASSLE; data munging is tricky.


2020.02.18
----------

Got the data formatted and ... it barfs. 


2020.08.18
----------

Have been on an analysis kick. Restricting data sets to those with N > 5 sites or more, and splitting the Joshua tree data by tree type, I can now get models to fit!

Here's the base model: geography + plant genetic distance explaining 

```
 Family: hurdle_gamma 
  Links: mu = log; shape = identity; hu = identity 
Formula: afst.x ~ pfst.x + geo.x + (1 | Pair) 
   Data: dat (Number of observations: 2328) 
Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
         total post-warmup samples = 4000

Group-Level Effects: 
~Pair (Number of levels: 17) 
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(Intercept)     1.26      0.26     0.88     1.91 1.01      706     1190

Population-Level Effects: 
          Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept    -2.18      0.32    -2.81    -1.58 1.01      517     1038
pfst.x        0.05      0.01     0.03     0.07 1.01     3131     2462
geo.x         0.32      0.04     0.24     0.40 1.00     3235     2700

Family Specific Parameters: 
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     2.23      0.06     2.11     2.36 1.00     3953     2714
hu        0.10      0.01     0.09     0.11 1.00     3649     2767

Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

Now, throw in interaction type (mutualism or antagonism):

```
 Family: hurdle_gamma 
  Links: mu = log; shape = identity; hu = identity 
Formula: afst.x ~ pfst.x + geo.x + Intxn.type + pfst.x:Intxn.type + (1 | Pair) 
   Data: filter(dat, Pair %in% keepers) (Number of observations: 2328) 
Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
         total post-warmup samples = 4000

Group-Level Effects: 
~Pair (Number of levels: 17) 
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(Intercept)     1.40      0.29     0.97     2.09 1.01     1010     1521

Population-Level Effects: 
                           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
Intercept                     -2.93      0.66    -4.26    -1.66 1.00     1588
pfst.x                         3.54      0.41     2.73     4.32 1.00     3050
geo.x                          0.30      0.04     0.22     0.38 1.00     4425
Intxn.typeMutualism            0.70      0.77    -0.81     2.22 1.00     1018
pfst.x:Intxn.typeMutualism    -3.50      0.41    -4.28    -2.68 1.00     3043
                           Tail_ESS
Intercept                      1938
pfst.x                         2471
geo.x                          2531
Intxn.typeMutualism            1633
pfst.x:Intxn.typeMutualism     2456

Family Specific Parameters: 
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     2.31      0.07     2.18     2.44 1.00     4870     2805
hu        0.10      0.01     0.09     0.11 1.00     5235     2597

Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

That's (1) a much bigger (nonzero) effect of plant genetics and (2) a nonzero interaction effect with mutualism, such that plant genetics has a smaller effect for mutualists?

```
 Family: hurdle_gamma 
  Links: mu = log; shape = identity; hu = identity 
Formula: afst.x ~ pfst.x + geo.x + Assoc.type + pfst.x:Assoc.type + (1 | Pair) 
   Data: dat (Number of observations: 2328) 
Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
         total post-warmup samples = 4000

Group-Level Effects: 
~Pair (Number of levels: 17) 
              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(Intercept)     1.29      0.26     0.89     1.91 1.00     1072     1694

Population-Level Effects: 
                         Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
Intercept                   -2.21      0.35    -2.93    -1.52 1.00     1019
pfst.x                       0.06      0.01     0.04     0.09 1.00     3600
geo.x                        0.36      0.04     0.27     0.44 1.00     4510
Assoc.typeMicrobe           -0.28      0.82    -1.89     1.36 1.00     1424
pfst.x:Assoc.typeMicrobe    -0.42      0.06    -0.52    -0.30 1.00     4076
                         Tail_ESS
Intercept                    1687
pfst.x                       2895
geo.x                        2969
Assoc.typeMicrobe            1726
pfst.x:Assoc.typeMicrobe     2874

Family Specific Parameters: 
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
shape     2.27      0.06     2.15     2.40 1.00     4601     2719
hu        0.10      0.01     0.09     0.11 1.00     4957     2564

Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

So *that* shows (1) much smaller effect of plant genetics, and (2) a significant interaction with microbe, such that the plant genetic effect is ... smaller? if the associate is a microbe. Hrm.


2020.08.20
----------

Model comparison among models (`.red` marks models "reduced" to remove an interaction):

```
Model comparisons:
              elpd_diff se_diff
mod.types       0.0       0.0  
mod.assoc     -21.1      19.5  
mod.types.red -41.0      17.3  
mod.assoc.red -41.5      17.8  
mod.base      -41.9      17.9  
```

That's solid preference for the model with interaction types.




