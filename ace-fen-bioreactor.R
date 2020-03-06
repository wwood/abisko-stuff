library(data.table)
library(ggplot2)
library(reshape2)
theme_set(theme_bw())
source("~/git/abisko-stuff/abisko.R")

add_bioreactor_metadata = function(abundance_in) {
    abundance = abundance_in
    abundance[, type := tstrsplit(sample,'\\.')[2]]
    abundance[, day1 := tstrsplit(sample,'\\.')[3]]
    abundance[, day := as.numeric(gsub('d','',day1))]
    abundance[, day1 := NULL]
    reboot_days = c(0, 28,63,105,147)
    reboot_round_names = factor(
        c('setup','setup','1','2','3'),
        levels=c('setup','1','2','3'))
    #abundance[, days_since_reboot := NULL]
    for (i in 1:length(reboot_days)) {
        abundance[day - reboot_days[i] >= 0, days_since_reboot := day-reboot_days[i]]
        abundance[day - reboot_days[i] >= 0, reboot_round := reboot_round_names[i]]
    }
    return(abundance)
}

add_duplicate_samples = function(abundance) {
    a3 = abundance[grep('duplicate',invert=T,sample)]
    added_c = a3[days_since_reboot==0 & reboot_round != 'setup' & type != 'feed']
    added_c[,type := 'baby_control']
    added_m = a3[days_since_reboot==0 & reboot_round != 'setup' & type != 'feed']
    added_m[,type := 'baby_monosaccharides']
    added_p = a3[days_since_reboot==0 & reboot_round != 'setup' & type != 'feed']
    added_p[,type := 'baby_polysaccharides']
    added = rbind(added_c, added_m, added_p)
    added[, sample := paste(sep='', sample, '.duplicate.',type)]
    return(rbind(a3, added))
}

ab = function(abundance) {
    return(add_duplicate_samples(add_bioreactor_metadata(abundance)))
}


read_community_profile = function() {
    d = fread('/srv/projects/abisko/aterrible_bins/18_a86_a62_promethion1_1529/20190424_epb1_flat_redo_vs_aterrible18_bins.fna.gz.coverm.csv')
    d2 = melt(d, id.vars=c('Genome')) # Complains about not all columns being the same - one is an int. Meh, annoying to see though.
    d2[, sample := gsub('.1.fq.gz .*','',gsub('.*/','',variable))]
    d2[, measure := gsub('.*.1.fq.gz ','',variable)]
    d2[, variable := NULL]
    d3 = dcast.data.table(d2, Genome+sample~measure, value.var='value')
    setnames(d3, c('Mean','Trimmed Mean'), c('mean_coverage','trimmed_mean'))
    setnames(d3, 'Relative Abundance (%)', 'relabund')
    setnames(d3, 'Covered Bases', 'covered_bases')
    setnames(d3, 'Genome','genome')

    ## a3 - incorporate duplicates so starting material for each is facetted, and mean_abundance
    a3 = add_duplicate_samples(add_bioreactor_metadata(d3))
    a3[, mean_abundance := mean(relabund), by='genome']

    ## Calculate growth rates
    a3[, genome_type_round := paste(genome,type,reboot_round)]

    ## Day 1
    current_round = 1
    prev_round = 0
    tm = merge(
        a3[days_since_reboot==prev_round,.(genome_type_round,relabund)],
        a3[days_since_reboot==current_round,.(genome_type_round,relabund)],
        by='genome_type_round'
    )[,.(gtr=genome_type_round, g=relabund.y - relabund.x)]
    a3[days_since_reboot==current_round, growth := tm[gtr==genome_type_round]$g, by=genome_type_round]
    ## Day 3
    current_round = 3
    prev_round = 1
    tm = merge(
        a3[days_since_reboot==prev_round,.(genome_type_round,relabund)],
        a3[days_since_reboot==current_round,.(genome_type_round,relabund)],
        by='genome_type_round'
    )[,.(gtr=genome_type_round, g=relabund.y - relabund.x)]
    a3[days_since_reboot==current_round, growth := tm[gtr==genome_type_round]$g, by=genome_type_round]
    ## Day 7
    current_round = 7
    prev_round = 3
    tm = merge(
        a3[days_since_reboot==prev_round,.(genome_type_round,relabund)],
        a3[days_since_reboot==current_round,.(genome_type_round,relabund)],
        by='genome_type_round'
    )[,.(gtr=genome_type_round, g=relabund.y - relabund.x)]
    a3[days_since_reboot==current_round, growth := tm[gtr==genome_type_round]$g, by=genome_type_round]
    ## Day 14
    current_round = 14
    prev_round = 7
    tm = merge(
        a3[days_since_reboot==prev_round,.(genome_type_round,relabund)],
        a3[days_since_reboot==current_round,.(genome_type_round,relabund)],
        by='genome_type_round'
    )[,.(gtr=genome_type_round, g=relabund.y - relabund.x)]
    a3[days_since_reboot==current_round, growth := tm[gtr==genome_type_round]$g, by=genome_type_round]

    a3[, type := factor(type, levels=c('parent','baby_control','baby_monosaccharides','baby_polysaccharides','feed'))]

    return(a3)
}

read_aterrible15_community_profile = function() {
    d = fread('bash -c \'cat <(head -1 /srv/projects/abisko/aterrible_bins/15_epb1_assembly86_assembly62/coverm_20190301_epb_flat.csv) <(tail -n+3 /srv/projects/abisko/aterrible_bins/15_epb1_assembly86_assembly62/coverm_20190301_epb_flat.csv)\'', header=T)
    setnames(d, 'Relative Abundance (%)', 'relative_abundance')
    d[, sample := gsub('.1.fq.gz$','',Sample)]
    d[, Sample := NULL]
    setnames(d, 'Genome', 'genome')
    setnames(d, c('Mean','Trimmed Mean'), c('mean_coverage','trimmed_mean'))

    d[grep('a62_',genome), new_genome_name := gsub('a62_','62_',genome)]
    d[grep('a86_',genome), new_genome_name := gsub('a86_','86_86_',genome)]

    return(d)
}

read_dereplicated_memberships = function () {
    d = fread('/srv/projects/abisko/aterrible_bins/18_a86_a62_promethion1_1529/clustering/dereplicated_memberships.csv')
    setnames(d, 'genome', 'genome_path')

    d[source=='old', genome := paste(sep='', 'old_',gsub('.fna$','',gsub('.*/','',genome_path)))]
    d[source=='62', genome := paste(sep='', '62_',gsub('.fna$','',gsub('.*/','',genome_path)))]
    d[source=='86', genome := paste(sep='', '86_',gsub('.fna$','',gsub('.*/','',genome_path)))]
    d[source=='mitch', genome := paste(sep='', 'mitch_',gsub('.fna$','',gsub('.*/','',genome_path)))]

    return(d)
}

read_coverm_strangeness_csvs = function() {
    files = c('/srv/projects/abisko/aterrible_bins/18_a86_a62_promethion1_1529/coverage_strangeness_investigation/alpha6_aterrible15.csv',
              '/srv/projects/abisko/aterrible_bins/18_a86_a62_promethion1_1529/coverage_strangeness_investigation/alpha4_aterrible15.csv',
              '/srv/projects/abisko/aterrible_bins/18_a86_a62_promethion1_1529/coverage_strangeness_investigation/alpha4_aterrible18.csv',
              '/srv/projects/abisko/aterrible_bins/18_a86_a62_promethion1_1529/coverage_strangeness_investigation/alpha6_aterrible18.csv')
    d = data.table(file=files)[,fread(file),by=file]
    d[, file := gsub('.*/','',file)]
    return(d)
}

read_gtdbtk = function() {
    a = fread('/srv/projects/abisko/aterrible_bins/18_a86_a62_promethion1_1529/gtdbtk4/gtdbtk.ar122.summary.tsv')
    b = fread('/srv/projects/abisko/aterrible_bins/18_a86_a62_promethion1_1529/gtdbtk4/gtdbtk.bac120.summary.tsv')
    return(
        rbind(a,b)[,.(genome=user_genome, taxonomy=classification)]
    )
}

read_idm = function() {
    d = fread('/srv/projects/abisko/shotgun_abundance/129_epb1_idm/idm_output/KO.matrix.tab',header=T)
    d2 = melt(d, id.vars='ID', variable.name='sample1', value.name='count')
    d2[, sample := gsub('.1.fq.gz.diamond.txt.gz','',sample1)]
    d2[, sample1 := NULL]
    setnames(d2, 'ID', 'ko')
    return(d2)
}

read_gc = function(dna_measured_days) {
    ## Notes here:
    ##
    ## Only the gas is measured for CH4 content, the dissolved is calculated from that using Henry's law
    ##
    ## Gas is sparged with N2 each time new peat is added to the parent reactor,
    ## but there might be dissolved gas still (not clear how much would come out
    ## of solution, but we assume all).
    ##
    ## CO2 has not yet been adjusted for pH of the liquid
    d = fread('cat \'/srv/projects/abisko/data/bioreactors/epb/operation data/Raw data csv/Parent GC.csv\'')
    m2 = t(d[,2:ncol(d),with=F])
    d2 = data.table(m2)
    setnames(d2, c('day_fraction','CH4_dissolved_produced','CO2_dissolved_produced',
                   'CH4_gas_produced','CO2_gas_produced'))
    d2[, day := round(day_fraction+0.5)] #=> always round up, because initial sample was late at night, i.e. after every other sample.
    d2[, type := 'parent']
    d2[, CH4_production := CH4_dissolved_produced + CH4_gas_produced]
    d2[, CO2_production := CO2_dissolved_produced + CO2_gas_produced]
    ## Checking
    stopifnot(!any(!(dna_measured_days %in% d2$day)))
    return(d2)
}

read_hplc_file = function(file, whitelist_chemicals, is_parent_reactor) {
    #print(file)
    #file = 'Parent HPLC.csv'
    d = fread(paste('/srv/projects/abisko/data/bioreactors/epb/operation data/Raw data csv/',file,sep=''),skip=1)

    d2 = d[3:nrow(d)]
    setnames(d2, as.character(c('chemical', t(d[1,2:ncol(d)]))))
    d2[, chemical := gsub(' .*','',chemical)]

    d3 = d2[chemical %in% whitelist_chemicals]
    stopifnot(length(unique(d3$chemical)) == length(whitelist_chemicals))

    feed_sample_ids = which(d[2]=='Feed') # Only in the Parent HPLC csv
    if (is_parent_reactor) {
        warning("The feed sample HPLC data has not been parsed from the parent HPLC CSV yet.")
        d3 = d3[,-feed_sample_ids,with=F]
    }

    d4 = data.table(melt(as.data.frame(d3), id.vars='chemical', variable.name='day_fraction', value.name='concentration', variable.factor=F))[!is.na(concentration)]
    d4[, day_fraction := as.numeric(as.character(day_fraction))]
    if (is_parent_reactor) {
        d4[, day := round(day_fraction+0.5)]
    } else {
        d4[, day := round(day_fraction)]
    }
    return(d4)
}

read_hplc = function(dna_measured_days) {
    whitelist_chemicals = c('Hexanoate','Valerate','nButyrate','Propionate','Acetate','Ethanol')
    d = data.table(
        file=c('Control HPLC.csv',
               'Mono HPLC.csv',
               'Parent HPLC.csv',
               'Poly HPLC.csv'),
        type=factor( # Factor must line up in all functions, expect R takes care if there's any issue tho
            c('baby_control','baby_monosaccharides','parent','baby_polysaccharides'),
            levels=c('parent','baby_control','baby_monosaccharides','baby_polysaccharides','feed')))

    d2 = d[, read_hplc_file(file, whitelist_chemicals, type=='parent'), by='file,type']
    d2[, file := NULL]

    ## Find the day_fraction cutoffs between before babies and after babies, use the cutoffs below
    ## d2[day==28 & type=='parent']
    ## d2[day==63 & type=='parent']
    ## d2[day==105 & type=='parent']
    ## d2[day==147 & type=='parent']
    d2[, reboot_round := 'setup']
    d2[day_fraction > 27.68, reboot_round := 'test']
    d2[day_fraction > 62.7, reboot_round := '1']
    d2[day_fraction > 104.6, reboot_round := '2']
    d2[day_fraction > 146.6, reboot_round := '3']

    ## Check if the days line up
    #data.table(table(d2[day %in% dna_measured_days]$day)) #=> Seems wrong there should be 24 or 12 for each day > 63 - caused by an error in the input spreadsheet - getting Rob to fix.
    warning("The baby hplc data is currently incomplete")

    d2[, sample := paste(sep='', 'r2.',type,'.d',day)]

    return(d2)
}

gather_dna_analysed_days = function() {
    dna_sample_names = fread('/srv/projects/abisko/aterrible_bins/18_a86_a62_promethion1_1529/samples_with_dna.list',header=F)
    dna_sample_names[, day := as.numeric(gsub('.*d','',V1))]
    return(dna_sample_names$day)
}

add_taxonomy_to_community_profile = function(community_profile, gtdbtk) {
    if (!is.null(community_profile$relabund)) {
        genome_abundance_order = unique(community_profile[,mean(relabund),by=genome][order(-V1)]$genome)
    } else {
        ## For use outside the local do.R e.g. in particles/mapping/15 we may not have relabund, so dummy something
        genome_abundance_order = unique(community_profile[order(genome)]$genome)
    }

    gtdbtk[domain=='d__', domain := '']
    gtdbtk[phylum=='p__', phylum := '']
    gtdbtk[class_name=='c__', class_name := '']
    gtdbtk[order_name=='o__', order_name := '']
    gtdbtk[family=='f__', family := '']
    gtdbtk[genus=='g__', genus := '']
    gtdbtk[species=='s__', species := '']
    gtdbtk[, tax := as.character(NA)]
    gtdbtk[species!='', tax := paste(genus, species)]
    ## Remove any label with a number in it as not a useful name
    gtdbtk[grep('[01-9]',tax), tax := NA]
    gtdbtk[is.na(tax) & genus!='', tax := genus]
    gtdbtk[grep('[01-9]',tax), tax := NA]
    gtdbtk[is.na(tax) & family!='', tax := family]
    gtdbtk[grep('[01-9]',tax), tax := NA]
    gtdbtk[is.na(tax) & order_name!='', tax := order_name]
    gtdbtk[grep('[01-9]',tax), tax := NA]
    gtdbtk[is.na(tax) & class_name!='', tax := class_name]
    gtdbtk[grep('[01-9]',tax), tax := NA]
    gtdbtk[is.na(tax) & phylum!='', tax := phylum]
    gtdbtk[grep('[01-9]',tax), tax := NA]
    gtdbtk[is.na(tax) & domain!='', tax := domain]
    gtdbtk[is.na(tax), tax := 'no domain assigned']
    gtdbtk[, class_tax := paste(class_name, tax)]
    gtdbtk[class_name == tax, class_tax := class_name]
    gtdbtk[, class_tax := gsub('c__','',class_tax)]

    community_profile = merge(community_profile, gtdbtk[,.(genome, class_tax)], by='genome', all.x=T)
    community_profile[, order_number := as.numeric(factor(genome, levels=genome_abundance_order))]
    community_profile[, name := paste(class_tax, order_number)]
    community_profile[genome == 'unmapped', name := 'unmapped']
    community_profile[, name := factor(name, levels=community_profile[order(order_number),unique(name)])]
    return(community_profile)
}


read_enrichm_ko_hmm = function() {
    d = melt(fread('/srv/projects/abisko/aterrible_bins/18_a86_a62_promethion1_1529/enrichm/enrichm_out_ko_hypothetical_cazy/ko_frequency_table.tsv'), id.vars='ID', variable.name='genome', value.name='ko_count')[ko_count > 0]
    setnames(d, 'ID', 'ko')
    return(d)
}

read_dna_extractions = function(types_factor_levels, community_profile) {
    d = fread('/srv/projects/abisko/data/bioreactors/epb/operation\ data/Extractions_Jamie-DNA.csv')[,1:11,with=F]
    ##d[Date %in% c('180903','181015')] #=> Not sure what these are - asking Rob.
    setnames(d, c('Day1','Month','Date','day','Type','Source','sample_order','Sample','extraction_order','dna_extraction_set','dna_yield'))
    d2 = d[,.(day,Type,sample_order,extraction_order,dna_extraction_set,dna_yield,rob_name=paste(Type,Date))]
    d2[,type := factor(NA, levels=types_factor_levels)]
    d2[Type == 'EPB', type := 'parent']
    d2[Type == 'Control', type := 'baby_control']
    d2[Type == 'Mono', type := 'baby_monosaccharides']
    d2[Type == 'Poly', type := 'baby_polysaccharides']
    d2[Type == 'Feed', type := 'feed']

    stopifnot(nrow(d2[is.na(type)])==2) # EPB_2 - waiting on Rob as above.
    d3 = d2[!is.na(type),.(day,type,sample_order,extraction_order,dna_extraction_set,dna_yield,rob_name)]
    d3[, sample := paste('r2.',type,'.d',day,sep='')]

    ## Check to make sure they line up - shouldn't have community profile from
    ## sample w/o extraction
    profiled_samples1 = community_profile[genome=='unmapped'][grep(invert=T,'duplicate',sample)]$sample
    stopifnot(length(profiled_samples1[!profiled_samples1 %in% d3$sample])==0)
    return(d3)
}


read_rna_extractions = function(types_factor_levels) {
    d = fread('/srv/projects/abisko/data/bioreactors/epb/operation\ data/Extractions_Jamie-RNA.csv')#[,1:11,with=F]
    setnames(d, c('Day1','Month','Date','day','CycleDay','Type','Source','sample_order',
                  'Unknown','Sample','Priority','Extracted','rna_yield','rin','rna_extraction_order',
                  'rna_extraction_notes'))
    d2 = d[,.(day,Type,sample_order,rna_extraction_order,rna_yield,rin)]
    d2[,type := factor(NA, levels=types_factor_levels)]
    d2[Type == 'EPB', type := 'parent']
    d2[Type == 'Control', type := 'baby_control']
    d2[Type == 'Mono', type := 'baby_monosaccharides']
    d2[Type == 'Poly', type := 'baby_polysaccharides']
    d2[Type == 'Feed', type := 'feed']

    stopifnot(nrow(d2[is.na(type)])==0)
    d3 = d2[!is.na(type),.(day,type,sample_order,rna_extraction_order,rna_yield,rin)]
    d3[, sample := paste('r2.',type,'.d',day,sep='')]

    return(d3[!is.na(rna_yield)])
}



#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

types_factor_levels = c('parent','baby_control','baby_monosaccharides','baby_polysaccharides','feed')
reboot_days = c(0, 28,63,105,147)

dna_measured_days = gather_dna_analysed_days()
community_profile1 = read_community_profile()
aterrible15_community_profile = read_aterrible15_community_profile()
dereplicated_memberships = read_dereplicated_memberships()
#coverm_strangeness_csvs = read_coverm_strangeness_csvs()
gtdbtk = split_taxonomy(read_gtdbtk(), has_root=F) #=> gtdbtk still running.
community_profile = add_taxonomy_to_community_profile(community_profile1, gtdbtk)
idm = read_idm()
gc = read_gc(dna_measured_days)
hplc = read_hplc(dna_measured_days)
enrichm_ko_hmm = read_enrichm_ko_hmm()
dna_extractions = read_dna_extractions(types_factor_levels, community_profile)
rna_extractions = read_rna_extractions(types_factor_levels)

## Waiting on Rob for
## 1. fix to hplc spreadsheet for babies
## 2. Gas production data for babies
## 3. HPLC data for final 2 points

