# Abisko-related code

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
