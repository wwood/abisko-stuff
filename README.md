**This code is specific for the [IsoGenie](https://isogenie.osu.edu/) project, and unlikely to be of more general use.**

### Abisko-related R code

To use in R, source it, something like so:
```
source('~/git/abisko-stuff/abisko.R')
```

Then there is functions to e.g. add metadata
```
$ R
> source('~/git/abisko-stuff/abisko.R')
> d = data.frame(sample=c('20120700_E1M','20110600_S1S','20110600_S1M'))
> d
        sample
1 20120700_E1M
2 20110600_S1S
3 20110600_S1M
> d2 = add_metadata(d)
> d2
        sample site year core month core_name depth site_ssplit
1 20120700_E1M    E 2012    1    07 201207_E1     M           E
2 20110600_S1S    S 2011    1    06 201106_S1     S   S_surface
3 20110600_S1M    S 2011    1    06 201106_S1     M S_anaerobic
```

And split taxonomy strings
```
$ R
> source('~/git/abisko-stuff/abisko.R')
> d = data.frame(taxonomy='Root;d__bacteria;p_phylum;came;over;for;great;spagbog')
> d2 = split_taxonomy(d)
> str(d2)
'data.frame':	1 obs. of  8 variables:
 $ taxonomy  : Factor w/ 1 level "Root;d__bacteria;p_phylum;came;over;for;great;spagbog": 1
 $ domain    : Factor w/ 1 level "d__bacteria": 1
 $ phylum    : Factor w/ 1 level "p_phylum": 1
 $ class_name: Factor w/ 1 level "came": 1
 $ order_name: Factor w/ 1 level "over": 1
 $ family    : Factor w/ 1 level "for": 1
 $ genus     : Factor w/ 1 level "great": 1
 $ species   : Factor w/ 1 level "spagbog": 1
```
If there is no `Root`, then use `split_taxonomy(d, has_root=F)`.
