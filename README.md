# Abisko-related code

To use in R, source it, something like so:
```
source('~/git/abisko-stuff/abisko.R')
```

Then there is functions to e.g. add metadata
```
$ R
> source('~/git/abisko-stuff/abisko.R')
> d = data.frame(sample=c('20120700_E1M','20110600_P1D'))
> d
        sample
1 20120700_E1M
2 20110600_P1D
> source('abisko.R')
> d2 = add_metadata(d)
> d2
        sample site year core month core_name depth site_ssplit
1 20120700_E1M    E 2012    1    07 201207_E1     M           E
2 20110600_P1D    P 2011    1    06 201106_P1     D           P
```
