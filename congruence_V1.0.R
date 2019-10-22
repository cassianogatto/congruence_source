# save this document as a text file in the wd folder as 'congruence_V1.0.R'
library(sf)
library(tidyverse)
library(ggplot2)
library(ggmap) # if using ggmap capabilities not presented in this tutorial)
# #_____________________________________________________________________________________#
# FUNCTION SP_LEVELS_BASED_ON_SIMILARITIES
'similarity levels: apply congruence threshold over a set of n fixed levels - CLOSED (if species are repeated in sucessive levels, or OPEN systems'
'plots all species at a threshold level to species focus in the first level round'
'the second level accounts for all species within min_sim accounted for all species in the first level'
'relies on function inter_patt to analyse and PLOT the maps'
'NEEDED split plot FALSE from save intersect to allow save intersect and union without plotting'
'from version 9.5 added map_species to inter_patt_levels'
'requires functions: inter_patt_level'
'Column names sim_df are species1, species2, C_thresh (for Csim)'
sp_levels_sim = function (sp_focus, sim_df = sim_df_ae,
                          map_species = mapAE_SA,
                          max_sim = 1, min_sim = 0.75, n_levels = 5,
                          save_intersect = TRUE,
                          new_window = FALSE, row_col = c(2,5), background_map = sa,
                          intersection_plot = TRUE,...) {
    # sp_focus is a species or a vector of species or a data frame with SCINAME column
    # print('sp_levels_sim_V1') # works with data frames (not only vectors); sp_df replaces spp_roll
    print('sp_levels_sim_V2.2') # intersection_plot = FALSE option  -> fast and clean closed/open system check
    # print("plots all species in a level of sim to species focus in the first level round; 
    #       the second level accounts for all species within min_sim accounted for all species in the first level, etc; 
    #       relies on function inter_patt to analyse and plot the maps; turn plot and feature calculus off with intersection_plot = FALSE")
    # V2.1 -> unified df with species and levels (species_level) + inter_union_ratio by level in a new slot
    # V2.2 -> updated and fixed intersection and union to be cumulative across levels - union was not working in levels>2
    # 
    # checks, data frames, create lists, etc
    if(is.vector(sp_focus)){  sp_df = data_frame('SCINAME'= sp_focus)  
    } else { if(is.data.frame(sp_focus)) { sp_df = sp_focus %>% select(SCINAME) } else {print('wrong type of sp input'); break  
    }}
    list_spp = list()
    unique_list_spp = data.frame()
    #check column Csim to C_thresh
    if ('Csim' %in% names(sim_df) & !('C_thresh' %in% names(sim_df)) ){sim_df = sim_df %>% select(species1,species2, C_thresh = Csim, everything())}
    
    list_spp[['info']] = paste(' sp_sim_levels V1 - congruence in a max number of levels (interactions)', '/n',
                               'pattern derived from', nrow(sp_df), 'species:', sp_focus, '/n',
                               'with max', max_sim, 'and min', min_sim, 'and ', n_levels, 'levels of analysis,','/n',
                               ifelse(intersection_plot,'with plot and calculus of intersections','with NO plot and calculus of intersections'))
    level_down = list_spp[['pattern']] # I didn't create the object before but it keeps just with the reference? hmmm how da hell is it workin'? Prob NOT! :/
    status = 'OPEN'
    if(intersection_plot){   if(new_window){ x11(); par(mar=c(2,5,2,2), mfrow = row_col)   }    }
    
    # LOOP LEVELS - N_LEVELS - remember the 'trick' : each time sp_df reload, it cames with
    #  the last list loaded...
    #  so its a feedback in itself and how the list grows in each cycle of the loop
    for (n in 1:n_levels){
        sp_df = sim_df %>% select(species1,species2,C_thresh) %>% filter(species1 %in% sp_df[,'SCINAME']) %>% 
            filter(as.numeric(C_thresh) >= as.numeric(min_sim) & as.numeric(C_thresh) <= max_sim) %>% 
            mutate('C_thresh' = round(as.numeric(C_thresh),3)) %>% arrange(desc(C_thresh)) %>%
            select(species1,species2) %>% gather %>% select(-key, 'SCINAME' = value) %>% unique %>% as.data.frame
        level = paste ('level', n, sep='')
        
        list_spp[[level]] = sp_df      # if someone wants individual level lists turn it on - it makes the final object a lit bit too confuse but makes things easier to plot
        
        print(paste(min_sim, ' level ', n))
        print(sp_df)
        # comparison between species in this with previous level - is the same set of species?
        if (all(unique(sp_df[,'SCINAME']) %in% unique(level_down[,'SCINAME']))) {
            status = 'CLOSED'
            list_spp[['status']] = paste ('CLOSED - after ', n-1,'levels',sep=' ')
            print(paste ('CLOSED - after ',n-1,'levels, with', length(unique(sp_df[,'SCINAME'])),'species',sep=' ')); break }
        
        list_spp[['species_level']] = rbind(list_spp[['species_level']], data.frame('SCINAME' = sp_df, 'level' = n)) %>% arrange(SCINAME,level)
        level_down = sp_df # backup copy for subsequent comparison
        title_fcn = paste (abbrev(sp_focus[1]),' ',round(min_sim,2), '\nlevel ', n, ', with ',
                           length(unique(sp_df[,'SCINAME'])),' species',sep='')
        #______________________________________________________________________________________________________________________________________
        # intersection & plots
        if (intersection_plot){     # if inter_patt = T  APPLY FUNCTION INTER_PATT inter_patt and PLOT maps and save the intersection calculated at this level
            
            inter_patt_level = inter_patt(sp_df = sp_df, map_species =  map_species, title = title_fcn, new_window = F, 
                                          background_map = background_map, plot_background = TRUE)
            
            # inter_patt is a list with 'intersect' & 'union'
            
            if (n > 1) { intersect_level = st_intersection(inter_patt_level[['intersect']] , intersect_level)  # make it cumulative across levels but I dont think is necessary - I already fixed the problem in inter_patt - pls check last version
            }  else { intersect_level = inter_patt_level[['intersect']] }
            
            assign(paste('inter_level', n, sep=''), intersect_level)    # modified 13Feb2019
            
            if (n > 1) { union_level = st_union(inter_patt_level[['union']], union_level)        # make it cumulative
            }  else  { union_level = inter_patt_level[['union']] }                    ## added 03/jan/19  # modified 12 fev 2019
            
            # calculus of intersection/union ratio per level # added 9/jan/2019
            list_spp[['inter_union_ratio']][[paste('level', n)]] = st_area(intersect_level)/ st_area(union_level)
            # calculus of the intersect area loss after new level species
            if (n > 1) {
                inter_last_level = get(paste('inter_level', n-1, sep=''))
                list_spp[['inter_loss']][[paste('level', n)]] = st_area(intersect_level)/ st_area(inter_last_level) }
            # break no intersection
            if (is.null(intersect_level)) { list_spp[['status']] = paste ('NO INTERSECTION on the level', n)
            print(paste ('NO INTERSECTION in the ', n,'level', sep=' '))
            break }
        } else {  print('No intersection evaluation and no plots. If needed pls use intersection_plot = TRUE')
            list_spp[['inter_loss']][[paste('level', n)]] = paste('no intersection calculus in this mode')        }
        #______________________________________________________________________________________________________________________________________    
        
        if (n == n_levels){ list_spp[['status']] = paste ('OPEN - after ',n,'levels with',length(unique(sp_df[,'SCINAME'])),'species',sep=' ')
        print(paste ('OPEN - after ', n,'levels with', length(unique(sp_df[,'SCINAME'])),'species', sep=' ')) }
        
        unique_list_spp =  list_spp[['species_level']] #rbind(unique_list_spp,sp_df) %>% unique
        
    } # for loop ends here
    # create a line to convert many lists to one data frame here
    # do.call("rbind", lapply(your_list, data.frame, stringsAsFactors=FALSE))
    if (intersection_plot){     if(save_intersect) { list_spp[['intersection']] = intersect_level; list_spp[['union']] = union_level }   }
    
    if(!intersection_plot){ list_spp = cbind(unique_list_spp, data.frame('min_sim'= rep(min_sim, nrow(unique_list_spp)),
                                                                         'n_lev'= rep(ifelse(status == 'CLOSED', n-1, n_levels),nrow(unique_list_spp)),
                                                                         'status'= rep(status, nrow(unique_list_spp)))) 
    list_spp = unique(list_spp) %>% arrange(desc(min_sim),level,SCINAME)
    }
    return(list_spp)
}
if(examples){
    test = sp_levels_sim(sp_focus = 'Capito auratus', sim_df = sim_df, max_sim = 1, min_sim = 0.6, n_levels = 4,
                         save_intersect = TRUE, new_window = TRUE, row_col = c(1,4), background_map = sa)
}

#___________________________________________________________________________________________________________________________#


#___________________________________________________________________________________________________________________________#
#
'find the congruence threshold between a closed and an open system for all species in a list'
'a data frame is generated with species in each level of congruence and status (open, closed), etc'
'use a data frame with species in a SCINAME column. tab_min excludes status==OPEN data...'
'requires functions: sp_levels_sim (requires inter_patt)'
'V9.5 removed "areas" and "tab_min" options - these transformations can be made outside via left_join with whatever area_df or summaries, etc.'

threshold_levels = function (spp_df, min_sim_init = .9, min_sim_step = .01, n_levels = 6,
                             sim_df = sim_df_lowlands, map_species = mapAE_SA,...){
    
    print('this is threshold_levels V1.3')# no areas or tab_min #print('this is threshold_levels V1.2')# V1.1 with minimal data frame (summarised) as result option #V 1.2 improvement in final df
    threshold_results = data.frame()
    threshold_tab_min = data.frame()
    min_sim_init_backup = min_sim_init
    if (!is.data.frame(spp_df)){print('spp_df needs to be a dataframe with species in SCINAME col'); break }
    if (all(names(spp_df) != 'SCINAME')){
        if (any(names(spp_df) == "species2")) { spp_df = spp_df %>% select('SCINAME' = species2)
        } else { spp_df = spp_df %>% select('SCINAME' = species1) }
    }
    
    for(i in 1: nrow(spp_df)){
        sp_focus = as.character(spp_df$SCINAME[i]);    sp_focus;      print(paste('progress ', i/nrow(spp_df)*100,'%'))
        #setup
        min_sim_init = min_sim_init_backup
        status = 'CLOSED'
        # loop status closed
        while(status == 'CLOSED'){
            if(min_sim_init < .1) {print('min_sim < .1'); break} # don't know why sometimes things go to negative loops minsim
            print(min_sim_init)
            
            #___call_function_sp_levels_sim _________________________________________________________________________________________
            list_spp_congruence_min = sp_levels_sim (sp_focus = sp_focus, sim_df = sim_df,
                                                     max_sim = 1, min_sim = min_sim_init, n_levels = n_levels,
                                                     intersection_plot = FALSE, map_species = map_species)
            #filter only the lower level for each congruence value
            # list_spp_congruence_min = list_spp_congruence_min %>% group_by(SCINAME) %>% filter(level == min(level)) %>% as.data.frame
            
            # while there is no spp > min_sim
            if (nrow(list_spp_congruence_min) == 0) { 
                min_sim_init = min_sim_init - min_sim_step; next  }  #criteria to decrease min_sim_init just works to intersection_plot = FALSE format
            
            #if OPEN break to next species without saving
            if (unique(list_spp_congruence_min$status)=='OPEN'){ break }
            if(nrow(list_spp_congruence_min) > 0){ list_spp_congruence_min = list_spp_congruence_min %>% filter(SCINAME != sp_focus) }
            
            # with species list > 0
            threshold_results_partial = cbind('sp_focus' = as.character(rep(sp_focus, nrow(list_spp_congruence_min))), list_spp_congruence_min)
            threshold_results_partial = threshold_results_partial %>% mutate_if(is.factor, as.character)
            
            threshold_results = rbind(threshold_results, threshold_results_partial)
            status = threshold_results_partial[nrow(threshold_results_partial), 'status'] 
            
            print('#'); print('##');  print('###')
            print(paste(sp_focus, min_sim_init, status)) # changed to status
            # status and min_sim_init set up for next round inner to while
            
            min_sim_init = min_sim_init - min_sim_step
            
            if(min_sim_init < .1) {print('min_sim < .1'); break}
            
            # > names(threshold_results)
            # [1] "sp_focus" "SCINAME"  "level"    "min_sim"  "n_lev"    "status" 
            
        } 
    }  # for loop ends here
    
    threshold_results = threshold_results %>% filter(status == 'CLOSED') %>%  filter(sp_focus != SCINAME) %>% 
        select(species1 = sp_focus, species2 = SCINAME, Cthres = min_sim, depth = level, status)
    
    return(threshold_results)
}

if(examples){
    l.serena = threshold_levels(spp_df = data.frame('species1'=c('Lepidothrix serena','Capito niger')), min_sim_step=.05 ,sim_df = sim_df_lowlands,tab_min=TRUE)
    
    thres_test = threshold_levels(spp_df = r.melano.73[1:3,], min_sim_init = .85, min_sim_step = .02, n_levels = 6)
}
#___________________________________________________________________________________________________________________________#

#___________________________________________________________________________________________________________________________#
'coherence function - plot one species with levels using plot_spp_inter; methods available'
'description:' 
{'This function takes a species as focus and plots species across a range
    of values of coherence; uses a similarity data frame, similarity max and min,
    and allows the choice of method (decimal (step=.1), user defined interval levels or fixed n levels).'
    'position = c("internal_to","overlap","contains") runs the internal function (sp_int_ext) and chooses
    only the selected categories - gives the status on the data frame printed and saved on object'
    'results in a list with species by level and if save_shapefiles saves final intersection and union'  
    'if  accum_intersection it will plot intersection and union accummulated across levels - otherwise 
    plot_spp_sim_maps will plot level inter and union'
    
    'dependencies: packages: tidiverse (dplyr), sf
    functions: plot_spp_sim_maps, bb_max_spp'}
coherence_to_sp = function(sp_focus, sim_df = sim_df_ae, max_sim = 1, min_sim = 0.2, map_species = mapAE_SA, 
                           method = 'decimal', #method = 'n_levels', n_levels = 5 ,  # method = 'fixed', fixed_levels = c(.99,.8,.7,.6), # default is the same as decimal
                           accum_intersection = TRUE,
                           position = c('overlap','contains', 'internal_to', NA), #  position = FALSE #(faster but with no internal external status)
                           n_levels = 5, # method = 'n_levels', # n_levels = trunc((max_sim-min_sim)*10), # option 
                           background_map = sa, row_col = c(2,4), 
                           save_shapefiles = TRUE, new_window = TRUE,
                           name_abbrev = TRUE, min_to_overlap = .95, iu_title = TRUE, ...  ){
    #print('coherence_to_sp V2 (internal function check')# print('coherence_to_sp V2.1 (internal function check') # modifications on the results[['species']] only one table
    #print('coherence_to_sp V2.2 (internal function check') # with Inters/Union ratio by level in data frame # print('coherence_to_sp V2.3') # acummulated_intersecion (and union) throughout levels (before was ploted only intra level)
    # print('coherence_to_sp V2.4') # background_plot deleted and modifications in the order of plot... trying to update IU to title in the same level, not lagged
    print('coherence_to_sp V2.5') # small adjusts in accum_intersection here and in plot_spp_sim_maps - vide comments above
    
    tab_spp = sim_df %>% filter(species1 == sp_focus) %>% filter (species2!= sp_focus) %>% 
        select(SCINAME = species2, C_thresh) %>% filter(C_thresh >= min_sim & C_thresh <= max_sim) %>%
        mutate_if(is.factor,as.character) %>% arrange(desc(C_thresh)) %>% as.data.frame
    if(nrow(tab_spp)==0) stop('unable to build a table with the current parameters')
    
    if (position[1] != FALSE){ # df internal external overlap (can control overlay criterion through 'min_to_overlap' if used at main function arguments)
        tab_spp1 = sp_int_ext (sp_focus = sp_focus, tab_spp = tab_spp, map = map_species,
                               min_to_overlap = ifelse(exists('min_to_overlap'),min_to_overlap,.95));
        tab_spp =  tab_spp1 %>% filter(status_sp_focus %in% position) %>% select(SCINAME, C_thresh, status_sp_focus, everything()) %>% as.data.frame() # !! everything
    }
    sp_map = map_species %>% filter(SCINAME == sp_focus)
    intersection = sp_map
    union = sp_map
    results = list()
    IU_ratio_accum_level = 1
    if(!exists('method')) { method = 'decimal' }
    results[['details']] = paste(sp_focus, min_sim,'-', max_sim,' method:', method, ifelse(method == 'n_levels', paste('n_levels=', n_levels),
                                                                                           ifelse(method == 'fixed', paste('levels=', fixed_levels), ifelse((method == 'decimal' | method == 'dec'), paste('.1 intervals'),''))))
    results[['species_list']] = data_frame()
    # print(results[['species_list']])
    bb_max_tab_spp = bb_max_spp(tab_spp, map_base = map_species)
    # print(bb_max_tab_spp)
    if(new_window){x11();par(mfrow = row_col, mar = c(3,3,3,2))}
    # define intervals of of sim by method and its parameters
    if (method == 'decimal' | method == 'dec') {  group_lvl = trunc(tab_spp$C_thresh * 10)  }
    if (method == 'n_levels' | method == 'n_lev') {    #group_lvl = split(tab_spp$SCINAME, cut( seq_along (tab_spp$SCINAME),n_levels, labels = FALSE)) # here just to consult options to split
        group_lvl = ntile(tab_spp$C_thresh, n_levels) } # same as above but simpler with ntile n_levels = 4 is 'quantile'
    if (method == 'fixed'){    group_lvl = cut( tab_spp$C_thresh,fixed_levels,labels=FALSE)  }
    ## plots with the function plot_spp_sim_maps
    for(ii in (as.numeric(unique(group_lvl)))){ # TEST tab_spp_lvl = tab_spp[which(group_lvl==ii),];print(tab_spp_lvl)}
        tab_spp_lvl = tab_spp[which(group_lvl==ii),]
        lvl_max_sim = tab_spp_lvl$C_thresh %>% max %>% round(.,3)
        lvl_min_sim = tab_spp_lvl$C_thresh %>% min %>% round(.,3)
        sp_abbrev = ifelse(name_abbrev, paste(str_sub(sp_focus,1,1),'.', word(sp_focus,2,2), sep = ''), sp_focus)
        # title & plot
        title = paste(sp_abbrev,'\n S:',lvl_max_sim,'-',lvl_min_sim,' N=', nrow(tab_spp_lvl),sep='') 
        print(title)
        print(tab_spp_lvl)
        
        plot(background_map$geometry, xlim = bb_max_tab_spp[c(1,3)], ylim = bb_max_tab_spp[c(2,4)])#, col=rgb(0,.1,0,1))
        
        ### Two colored title :)
        if(iu_title){
            if(exists('IU_ratio_accum_level')) {
                IU_title = paste('last_IU_ratio = ',round(as.numeric(IU_ratio_accum_level),2), sep='')
                if(as.numeric(IU_ratio_accum_level) < .5) { 
                    title(paste(title,' ',IU_title), col.main = 'dark red')
                } else {
                    title(paste(title,' ',IU_title), col.main = 'dark blue') 
                }}} else {title(title)}
        
        ## select level species and apply function plot_sim_spp_maps
        Pv_plot = plot_spp_sim_maps(tab_spp = tab_spp_lvl, background_plot = FALSE, map_species = map_species, 
                                    plot_level_inter_union = !accum_intersection)
        
        # if (st_crs(intersection) != st_crs(map_species)) {st_crs(intersection) = st_crs(map_species); st_crs(union) = st_crs(map_species) }
        # if (st_crs(Pv_plot[['intersection']]) != st_crs(map_species)) {st_crs(Pv_plot[['intersection']]) = st_crs(map_species) }
        # if (st_crs(Pv_plot[['union']]) != st_crs(map_species)){ st_crs(Pv_plot[['union']]) = st_crs(map_species)}
        
        if(IU_ratio_accum_level != 0) {intersection = intersection %>% st_intersection (Pv_plot[['intersection']])}
        union = union %>% st_union(Pv_plot[['union']])
        if (length(st_area(intersection)) == 0 ){IU_ratio_accum_level = 0 } else {
            IU_ratio_accum_level = as.numeric(st_area(intersection))/ as.numeric(st_area(union)) %>% as.numeric }
        
        plot(sp_map %>% st_geometry, lwd=1, add=T, col = rgb(0,0,0,.2))
        
        #results[['species_list']][[paste('LEVEL_',ii,' S:', lvl_max_sim,'-', lvl_min_sim, sep='')]] = tab_spp_lvl
        data = cbind(tab_spp_lvl, data_frame ('max_sim_level' = rep(lvl_max_sim, nrow(tab_spp_lvl)),
                                              'min_sim_level' = rep(lvl_min_sim, nrow(tab_spp_lvl)),
                                              'IU_accum_level' = rep(as.numeric(IU_ratio_accum_level), nrow(tab_spp_lvl))))
        results[['species_list']] = rbind(results[['species_list']],data)
        # saving only level values, NOT accumulated ones... needs to apply something like 'lapply function(x) inter = inter st_inter(x)' etc...
        if(save_shapefiles){results[['shapefiles']][[paste('level',ii)]] = list('intersection' = Pv_plot[['intersection']],'union' = Pv_plot[['union']])  }
        
        # plot intersection and union accumulated
        if(accum_intersection){ 
            plot(st_geometry(union), add = T, col = rgb(0,1,1,.01), lwd = 2, lty = 1)
            if(length(st_area(intersection)) != 0) {plot(st_geometry(intersection), add = T, col = rgb(.5,.5,0,.7), lwd = 2)}
        }
    }
    if(save_shapefiles){results[['intersection']] = intersection; results[['union']] = union}
    if(position[1] != FALSE) {results[['species_list']] = results[['species_list']] %>% select(sp_focus, status_sp_focus, SCINAME, 
                                                                                               C_thresh, IU_sp_sp, area_sp_sp_ratio, IU_accum_level, max_sim_level, min_sim_level, everything() )}
    results
}
# examples  coherence (need to updata to accum inter union at level or cummulative...etc just first and second are)
if(examples){
    r.melanos_coher_accum = coherence_to_sp(sp_focus = 'Rhegmatorhina melanosticta',sim_df = sim_df, map_species= mapAE_SA, max_sim = 1, min_sim = .5,
                                            method = 'n_levels', n_levels = 6, row_col = c(1,1), accum_intersection = TRUE, position = c('contains','internal_to','overlap',NA),
                                            background_map = sa, save_shapefiles = FALSE, new_window = TRUE,
                                            name_abbrev = TRUE, min_to_overlap = .95)
    r.melanos_coher_level = coherence_to_sp(sp_focus = 'Rhegmatorhina melanosticta',sim_df = sim_df, map_species= mapAE_SA, max_sim = 1, min_sim = .6,
                                            method = 'n_levels', n_levels = 6, row_col = c(2,3), accum_intersection = FALSE, position = c('contains','internal_to','overlap',NA),
                                            background_map = sa, save_shapefiles = FALSE, new_window = TRUE,
                                            name_abbrev = TRUE, min_to_overlap = .95)
    
    # TEST DECIMAL .9,.8,.7 etc
    coherence_to_sp(sp_focus = 'Rhegmatorhina melanosticta',min_sim = .6,max_sim = 1, method = 'dec',
                    sim_df = sim_df,  background_map = sa, map_species = mapAE_SA, new_window = T)
    # TEST fixed vector - interval levels
    coherence_to_sp(sp_focus = 'Xiphorhynchus pardalotus',min_sim = .6,max_sim = 1, method = 'fixed',
                    fixed_levels = c(.99,.85,.7,.65,.6),
                    sim_df = sim_df,  background_map = sa, map_species = mapAE_SA)
    # TEST n levels
    coherence_to_sp(sp_focus = 'Phlegopsis nigromaculata',min_sim = .5,max_sim = 1, method = 'n_levels',
                    n_levels = 8, row_col_par = c(2,4),
                    sim_df = sim_df, background_map = sa, map_species = mapAE_SA)
    # TEST with position
    coherence_to_sp(sp_focus = 'Rhegmatorhina melanosticta',min_sim = .6,max_sim = 1, method = 'dec',
                    sim_df = sim_df,  background_map = sa, map_species = mapAE_SA, new_window = T)
    # TEST with save_shapefiles = T
    test=coherence_to_sp(sp_focus = 'Capito auratus',min_sim = .5,max_sim = 1, method = 'n_levels', n_levels=5,
                         sim_df = sim_df, position = c('overlap','contains'),#, 'internal_to', NA),
                         background_map = sa, map_species = mapAE_SA, new_window = T,accum_intersection = TRUE,
                         save_shapefiles = TRUE)
    
}
#___________________________________________________________________________________________________________________________#
# sp_int_ext V1 - internal external
'to calculate internal or external status of sp_focus compared to other polygons in a df$SCINAME:'
'description'
{'overlap -> when the congruence is >=.95 (default) or any value min_to_overlap'
    'contains -> when the species focus contains > .98 of sp2 range'
    'internal_to -> sp1 is >=.98 inside of sp2 range (overlay = 0.98 of sp1 range)'
    'NA -> none of the above criteria applies'
    "to filter use (.) %>% filter(status_sp_focus == 'internal_to' | status_sp_focus == 'contains') %>% as.data.frame()"
    'V.9.5 map to map_species'}
sp_int_ext = function(sp_focus, tab_spp, map_species = mapAE_SA, 
                      min_to_overlap = .9, min_to_internal_to = .95, min_to_contains = .95,  # <<< --- parameters to internal, contains, overlap and NA !!!
                      full_df = FALSE,...){
    # print('sp_int_ext V1')
    # print('sp_int_ext V1.1') # min to contains and internal  - all defaults to .9
    print('sp_int_ext V1.2') # IU ratio sp by sp / ratio sp by sp areas
    if(is.vector(tab_spp)) { tab_spp = data_frame('SCINAME' = tab_spp) }
    sp_map = map_species %>% filter(SCINAME == sp_focus)
    area_sp_map = st_area(sp_map)
    df_results = data_frame()
    for (i in 1:length(tab_spp$SCINAME)){
        map_sp2 = map_species %>% filter(SCINAME == tab_spp$SCINAME[i])
        area_sp2 = st_area(map_sp2)
        inter_area = st_intersection(sp_map,map_sp2) %>% st_area
        ov_1 = inter_area/area_sp_map
        ov_2 = inter_area/area_sp2
        area_union = st_area(st_geometry(st_union(sp_map,map_sp2)))
        # status verification
        if (as.numeric(ov_1*ov_2) >= min_to_overlap) {status_sp1 = 'overlap'} else {
            if (as.numeric(ov_1) > as.numeric(ov_2)){
                if (as.numeric(ov_1) >= min_to_internal_to) {status_sp1 = 'internal_to'} else {status_sp1 = NA}} else {
                    if (as.numeric(ov_2) > min_to_contains) {status_sp1 = 'contains'}  else {status_sp1 = NA}}
        }
        # write df
        if(full_df){df_results = rbind(df_results,data_frame('sp_focus' = sp_focus, 'SCINAME' = tab_spp$SCINAME[i],
                                                             'status_sp_focus'= status_sp1, 'overlay' = inter_area,
                                                             'C_thresh' = as.numeric(round(ov_1*ov_2,3)),
                                                             'IU_sp_sp' = as.numeric(round(inter_area/area_union,3)),
                                                             'area_sp_sp_ratio' = as.numeric(round(area_sp_map/area_sp2,3)),
                                                             'area_sp1' = area_sp_map, 'area_sp2' = area_sp2  
        ))
        } else {df_results = rbind(df_results,data_frame('sp_focus' = sp_focus,'status_sp_focus'= status_sp1,
                                                         'SCINAME' = tab_spp$SCINAME[i],'C_thresh' = as.numeric(round(ov_1*ov_2,4)),
                                                         'IU_sp_sp' = as.numeric(round(inter_area/area_union,3)),
                                                         'area_sp_sp_ratio' = as.numeric(round(area_sp_map/area_sp2,3))
        )) } 
    }
    df_results %>% as.data.frame
}


# #___________________________________________________________________________________________________________________________#
# #     CO-OCCURRENCE PATTERNS FUNCTION  
# co_occurrence_patterns = function (test_mode = FALSE, D_df=D_df,   # similarity matrix (vide C_thresh function)
#                                    min_sim=0.7, # minimal similarity in first level choice
#                                    min_sim_vect = 0.7,  # figure out the exact meaning of this parameter!
#                                    min_sim_to_pattern = 0.5, # min_sim to first species
#                                    ratio_IU = 0.5,
#                                    method = 'default',#vector_sim', # methods = 
#                                    vector_sim=vector_sim,
#                                    plot_map=TRUE,
#                                    pdf_plot = FALSE,
#                                    plot_text = FALSE,
#                                    map = mapAE_SA, # map source to matrix species (same nammes)
#                                    save_spat_obj = FALSE, 
#                                    background = T,
#                                    background_map = neotropical,
#                                    spatial_check = TRUE,
#                                    valid_chk = TRUE,
#                                    abbreviated = FALSE, #excl
#                                    dir_shp = "C:/SIG2018/CO_OCCUR_PATT/patterns_shp",
#                                    dir_output= "C:/SIG2018/CO_OCCUR_PATT/patterns_df/"
# ) {
#     
#     library(dplyr)
#     library(sf)
#     oldw <- getOption("warn")
#     options(warn = -1)
#     
#     if (spatial_check) source('check_spatial_function.R') # loads check spat funct to check_spat
#     if (abbreviated) source('make_cep_names_function.R')  # cep names if needed
#     
#     # TEST MODE ONLY
#     if (test_mode) {D_df= sim_matrix[1000:2000,1000:2000]# default + test parameters
#     min_sim=0.7; method = 'default'; map=mapAE_SA; plot_map=TRUE; plot_text = FALSE; ratio_IU = 0.5;
#     min_sim_vect = 0.7; min_sim_to_pattern = 0.5; save_spat_obj = FALSE; abbreviated = FALSE;
#     background = TRUE; background_map = neo_merge; spatial_check = TRUE;  valid_chk = TRUE;  pdf_plot = TRUE;
#     dir_shp = "C:/SIG2018/CO_OCCUR_PATT/teste_jun2018";
#     dir_output= "C:/SIG2018/CO_OCCUR_PATT/teste_jun2018/"
#     }
#     output = list() # output
#     output[['meta_info']] = paste('Species co-occurrence patterns based on geographical distribution similarity', '\n', 'Presets:','\n',
#                                   'minimal similarity = ', min_sim, '\n','min similarity to build vector of similarity between species = ',min_sim_vect, '\n',
#                                   'min similarity to first(pattern) species = ', min_sim_to_pattern, '\n','min ratio between cumulative intersection and union = ', ratio_IU, '\n'  ) # meta-info
#     # [1] Matrix setup, markers, etc.
#     if (any(colnames(D_df) != rownames(D_df))){stop("Error: matrix not symmetrical|rows and columns are not identical")}
#     # Dist similarity table manipulation & filtering
#     for(z in 1:nrow(D_df)){D_df[z,z] = NA} # transform all 1 in NA for equal species in row/col
#     D_df = ifelse(D_df < min_sim, NA, D_df)   # all D_df < min_sim are converted to NA
#     filt_r_min = apply(D_df,1,function(x) all(is.na(x)))  # all values (spp similarities) < min_sim? or NA?
#     D_df = D_df[!filt_r_min,!filt_r_min]   # only spp with at least one similarity > min
#     if (any(colnames(D_df) != rownames(D_df))){stop("Error: matrix
#                                                     not symmetrical after filtering! Check filters to NA and min_sim")}
#     nome = rownames(D_df)
#     ratio_inter_union = 1
#     if(abbreviated){tab_abrev = data.frame('species' = as.character(map[["SCINAME"]]),  ########## species column identifyer
#                                            'sp_abrev' =  make.cepnames1(droplevels(map$SCINAME),nchar=6))} else { tab_abrev = data.frame('species' = as.character(map[["SCINAME"]])) }
#     if (test_mode) head(tab_abrev)
#     # vector_D_sim
#     if (method == 'vector_sim'){vector_D_sim = vector_sim} else {
#         # ordered vector with D_sim SUM (or means) by row
#         # vector_D_sim = order(rowMeans(D_df, na.rm = TRUE), decreasing = TRUE) # means or sums of >min values
#         vector_D_sim = order(rowSums(D_df, na.rm = TRUE), decreasing = TRUE)  ; if (test_mode) vector_D_sim
#     }
#     ###__ plot, graphical parameters  #
#     if (plot_map) {  x11()
#         #par(mfrow = c(1,2))#, mar = c(2,4,2,4)) 
#     }
#     
#     ### [2] MAIN LOOP - PATTERN LOOP _____###
#     ###___________________________________________________________________________________
#     
#     for (sp in vector_D_sim) {           # if(test_mode) sp=vector_D_sim ['X'] 
#         
#         # sorted vector of non NA's similarities in decreasing order
#         vector_Dsp = sort(D_df[sp,][!is.na(D_df[sp,])],decreasing = T)  #ORDERED VECTOR
#         paste('sp = ',sp, ' / ', rownames(D_df)[sp], sep = '')
#         
#         ####   first species to compare with sp == max similarity
#         sp_2 = which (D_df[sp,] == max(D_df[sp,],na.rm=T))[1] # [1] required if there are more than one spp with max D_sim
#         paste('sp_2 = ',sp_2, ' / ', rownames(D_df)[sp_2], sep = '')
#         
#         # DISTANCE SIMILARITY SP x SP_2
#         D_loop = D_df[sp,sp_2]; if (test_mode) D_loop
#         D_init = D_loop
#         
#         # Pattern configuration
#         sp_pattern = sp # keeps the #1 pattern species index
#         pattern_out = paste('pattern_',nome[sp_pattern],sep="")
#         output[[pattern_out]] = list()
#         
#         #loop parameters
#         loop = 'go'
#         n_loop = 0
#         nxt = 0
#         ratio_inter_union = 1
#         cat('                      PATTERN                      ')
#         cat('               ',pattern_out," Starts               ! ! 8==> O-:} ")
#         check_spat_results = list()
#         
#         # creation of empty inter and union pattern spatial objects
#         inter_comb = st_sf(id = 1, geometry = st_sfc(lapply(1, function(x) st_multipolygon()))); if (test_mode) inter_comb
#         union_comb = st_sf(id = 1, geometry = st_sfc(lapply(1, function(x) st_multipolygon()))); if (test_mode) union_comb
#         
#         # saving pattern plots pdf_plot = TRUE
#         # n_dev_copy = n_dev_copy + 1
#         
#         ###______ INSIDE PATTERN  - INNER LOOP -> WHILE / REPEAT ______###
#         #____________________________________________________________________________________________
#         repeat {
#             #____________________________________________________________________________________________      
#             # Setup vector D_sp / D_sp_2
#             n_loop = n_loop + 1
#             # ORDERED VECTOR 
#             vector_Dsp_2= sort(D_df[sp_2,][!is.na(D_df[sp_2,])],decreasing = T) 
#             if (test_mode) vector_Dsp_2 
#             
#             # if sp changes after first loop, recalculates vector_Dsp
#             if (sp != sp_pattern) { #sp_col = which(!is.na(D_df[sp,]), arr.ind = T)  
#                 vector_Dsp = sort(D_df[sp,][!is.na(D_df[sp,])],decreasing = T) #ORDERED VECTOR
#             }
#             if (test_mode) vector_Dsp
#             
#             ###______ Data Frame  Update  ________## transform to a function!
#             # out_table = 'unchanged'
#             if (n_loop == 1) { #fist loop of the pattern ? -> make the complete table...
#                 spec = c(rep(nome[sp],length(vector_Dsp)),rep(nome[sp_2],length(vector_Dsp_2)))
#                 spec_ind = c(rep(sp,length(vector_Dsp)),rep(sp_2,length(vector_Dsp_2)))
#                 spec_2 = c(names(vector_Dsp),names(vector_Dsp_2))
#                 spec_ind_2 = c( sapply(names(vector_Dsp),function(x)which(rownames(D_df) == x)),
#                                 sapply(names(vector_Dsp_2),function(x)which(rownames(D_df) == x)))
#                 D_sim = as.numeric(c(vector_Dsp, vector_Dsp_2))
#                 used = numeric(length(vector_Dsp)+length(vector_Dsp_2))
#                 
#             } else {  
#                 #  compare sp and sp_2 only to first column'spec'  
#                 if (!(sp %in% output[[pattern_out]][['df']][,'spec_ind']) |
#                     !(sp_2 %in% output[[pattern_out]][['df']][,'spec_ind'])) { # sp OR sp_2 %in% output[['df']][['spec_ind']]
#                     # sp or sp_2   # sp in spec but sp_2 not
#                     if ((sp %in% output[[pattern_out]][['df']][,'spec_ind']) &
#                         !(sp_2 %in% output[[pattern_out]][['df']][,'spec_ind'])) {  # sp in but sp_2 not in
#                         spec = c(rep(nome[sp_2],length(vector_Dsp_2)))
#                         spec_ind = c(rep(sp_2,length(vector_Dsp_2)))
#                         spec_2 = names(vector_Dsp_2)
#                         spec_ind_2 = sapply(names(vector_Dsp_2),function(x)which(rownames(D_df) == x))
#                         D_sim = as.numeric(vector_Dsp_2)
#                         used = numeric(length = length(vector_Dsp_2))
#                     }
#                     # sp_2 in spec but not sp
#                     if (!(sp %in% output[[pattern_out]][['df']][,'spec_ind']) &
#                         (sp_2 %in% output[[pattern_out]][['df']][,'spec_ind'])) {  #sp_2 in but sp not in
#                         spec = c(rep(nome[sp],length(vector_Dsp)))
#                         spec_ind = c(rep(sp,length(vector_Dsp)))
#                         spec_2 = names(vector_Dsp)
#                         spec_ind_2 = sapply(names(vector_Dsp),function(x)which(rownames(D_df) == x))
#                         D_sim = as.numeric(vector_Dsp)
#                         used = numeric(length = length(vector_Dsp))
#                     }
#                     # sp & sp_2 NOT in output same as n_loop == 1
#                     if (!(sp %in% output[[pattern_out]][['df']][,'spec_ind']) &
#                         !(sp_2 %in% output[[pattern_out]][['df']][,'spec_ind'])) { 
#                         spec = c(rep(nome[sp],length(vector_Dsp)),rep(nome[sp_2],length(vector_Dsp_2)))
#                         spec_ind = c(rep(sp,length(vector_Dsp)),rep(sp_2,length(vector_Dsp_2)))
#                         spec_2 = c(names(vector_Dsp),names(vector_Dsp_2))
#                         spec_ind_2 = c( sapply(names(vector_Dsp),function(x)which(rownames(D_df) == x)),
#                                         sapply(names(vector_Dsp_2),function(x)which(rownames(D_df) == x)))
#                         D_sim = as.numeric(c(vector_Dsp, vector_Dsp_2))
#                         used = numeric(length = length(vector_Dsp)+length(length(vector_Dsp_2)))
#                     }
#                 } # sp OR sp_2 %in% 
#             }   
#             tab = data.frame(spec,spec_2,spec_ind,spec_ind_2,D_sim = as.numeric(D_sim),used=as.numeric(used))
#             rownames(tab)= NULL; if (test_mode)tab 
#             output[[pattern_out]][['df']] = unique(rbind(as.data.frame
#                                                          (output[[pattern_out]][['df']]),tab)) 
#             
#             ###________ MAPS CALCULUS_POLYGONS 
#             # map sp sp_2 sp_pattern
#             if (n_loop ==1) { if (!abbreviated){ 
#                 ap = which(tab_abrev$species == nome[sp_pattern])} else {
#                     ap = which(tab_abrev$sp_abrev == nome[sp_pattern])    }   
#                 map_sp_pattern = map[ap,]
#                 map_sp_pattern = map_sp_pattern %>% group_by(SCINAME) %>% summarise() }
#             if (test_mode) map_sp_pattern
#             ab = ifelse(abbreviated,which(tab_abrev$sp_abrev == nome[sp]),
#                         which(tab_abrev$species == nome[sp]))
#             map_sp = map[ab,]
#             map_sp = map_sp %>% group_by(SCINAME) %>% summarise()# %>% st_cast("POLYGON")
#             map_sp = map_sp %>% st_buffer(0)
#             if (test_mode)map_sp
#             
#             ab_2 = ifelse(abbreviated,which(tab_abrev$sp_abrev == nome[sp_2]),
#                           which(tab_abrev$species == nome[sp_2]))
#             map_sp_2 = map[ab_2,] #%>% st_cast("MULTIPOLYGON")
#             map_sp_2 = map_sp_2 %>% group_by(SCINAME) %>% summarise()
#             map_sp_2 = map_sp_2 %>% st_buffer(0)
#             if (test_mode) map_sp_2
#             #check st_is_valid and fix with buffer...
#             if(valid_chk){
#                 if ((!st_is_valid(map_sp_pattern)) | (!st_is_valid(map_sp)) | (!st_is_valid(map_sp_2))){
#                     if (!st_is_valid(map_sp_pattern)) {map_sp_pattern = map_sp %>% st_buffer(0.05) %>% st_buffer(-0.05)}
#                     if (!st_is_valid(map_sp)) {map_sp = map_sp %>% st_buffer(0.05) %>% st_buffer(-0.05)}
#                     if (!st_is_valid(map_sp_2)) {map_sp_2 = map_sp_2 %>% st_buffer(0.05) %>% st_buffer(-0.05)}
#                 } else{ if (test_mode) cat('valid_chk shp sp, sp_2, sp_pattern => OK')
#                 }}
#             ### check st_geometry_type !=MULTIPOLYGON and force with st_cast...      # I'll try to run without these fixing stuff!
#             # if ((st_geometry_type(map_sp_pattern) != "MULTIPOLYGON") |
#             #     (st_geometry_type(map_sp) != "MULTIPOLYGON") |
#             #     (st_geometry_type(map_sp_2) != "MULTIPOLYGON")){
#             #   if (st_geometry_type(map_sp_pattern) != "MULTIPOLYGON") {
#             #     map_sp_pattern = map_sp_pattern %>%  st_cast("MULTIPOLYGON")}
#             #   if (st_geometry_type(map_sp) != "MULTIPOLYGON") {
#             #     map_sp = map_sp  %>% st_cast("MULTIPOLYGON")}
#             #   if (st_geometry_type(map_sp_2) != "MULTIPOLYGON") {
#             #     map_sp_2 = map_sp_2  %>% st_cast("MULTIPOLYGON")     }      }
#             
#             # EXCLUDED BY NOW
#             # Intersection & union of the Pattern
#             inter_sp = st_intersection(st_geometry(map_sp), st_geometry(map_sp_2)) %>% st_buffer(0)
#             union_sp = st_union(st_geometry(map_sp), st_geometry(map_sp_2)) %>% st_buffer(0)
#             if (test_mode)inter_sp
#             if (test_mode)union_sp
#             # check for erros and fixing forcing to multipolygon to avoid errors like: 'Error in as(st_geometry(x), "Spatial") : 
#             if (valid_chk){
#                 if ((!st_is_valid(inter_sp))|(!st_is_valid(union_sp))){ # forcing to multipolygon to avoid errors like: 'Error in as(st_geometry(x), "Spatial") : 
#                     if (!st_is_valid(inter_sp)){ inter_sp = st_buffer(inter_sp,0) }
#                     if (!st_is_valid(union_sp)){ union_sp = st_buffer(union_sp,0) }
#                 } else { if (test_mode) cat ('valid_chk shp inter & union -> OK')
#                 }}
#             if ((st_geometry_type(inter_sp) != "MULTIPOLYGON") |  (st_geometry_type(union_sp) != "MULTIPOLYGON")) {
#                 if (st_geometry_type(inter_sp) != "MULTIPOLYGON") {inter_sp = inter_sp %>% st_cast("MULTIPOLYGON") }
#                 if (st_geometry_type(union_sp) != "MULTIPOLYGON") {union_sp = union_sp %>% st_cast("MULTIPOLYGON") }
#             }
#             # CRS
#             st_crs(union_comb) = st_crs(inter_comb) = st_crs(map_sp)
#             # COMBINE: inter and union comb
#             if (n_loop == 1) {inter_comb = inter_sp; union_comb = union_sp
#             } else {
#                 inter_comb = st_intersection(st_geometry(inter_comb), st_geometry(inter_sp)) %>% st_buffer(0)
#                 union_comb = st_union(st_geometry(union_comb), st_geometry(union_sp)) %>% st_buffer(0)
#                 # Intersection & union of the Pattern
#                 # check for erros and fixing forcing to multipolygon to avoid errors like: 'Error in as(st_geometry(x), "combatial") : 
#                 if (valid_chk){
#                     if ((!st_is_valid(inter_comb))|(!st_is_valid(union_comb))){ # forcing to multipolygon to avoid errors like: 'Error in as(st_geometry(x), "combatial") : 
#                         # checkin and fixing polygon validity N2
#                         if (!st_is_valid(inter_comb)) {inter_comb = st_buffer(inter_comb,0)}
#                         if (!st_is_valid(union_comb)) {union_comb = st_buffer(union_comb,0)}
#                     } else { if (test_mode) cat('valid_chk shp inter & union combined are OK')
#                     }
#                 }
#                 if ((st_geometry_type(inter_comb) != "MULTIPOLYGON") | 
#                     (st_geometry_type(union_comb) != "MULTIPOLYGON")) {
#                     if (st_geometry_type(inter_comb) != "MULTIPOLYGON") {inter_comb = inter_comb %>% st_cast("MULTIPOLYGON") }
#                     if (st_geometry_type(union_comb) != "MULTIPOLYGON") {union_comb = union_comb %>% st_cast("MULTIPOLYGON") }
#                 }
#             }
#             if (test_mode)inter_comb
#             if (test_mode)union_comb
#             # ratio used to break pattern
#             ratio_inter_union = as.numeric(as.character(st_area(inter_comb)/st_area(union_comb)))
#             if (test_mode)ratio_inter_union
#             
#             ##_________CHECK LOOP _____________##
#             # condition on ratio inter/union
#             # leave the loop condition
#             if (ratio_inter_union < min_sim) {break; cat('Calculus polygons: ratio_inter_union < min_sim')}
#             
#             ####### first check spatial
#             if (spatial_check){check_spat_results = check_spatial()
#             if(check_spat_results[['status']]==FALSE){break; print(unlist(check_spat_results))}
#             if (test_mode) check_spat_results[['status']]
#             }
#             
#             ##___________ PLOTS ___________________##
#             #
#             if (plot_map){
#                 if(n_loop ==1){ 
#                     if (background){ 
#                         plot(union_sp,col=rgb(0,0,1,.2),  add=F,main=paste(nome[sp],' x ', nome[sp_2]))
#                         plot(st_geometry(background_map),lwd=1.5,add=T) #main=paste(gsub('pattern_','',pattern_out))
#                         
#                     } else {
#                         plot(union_sp,col=rgb(0,0,1,.2), xlim = x, ylim = y,
#                              main= paste(gsub('pattern_','',pattern_out)))
#                     }
#                     plot(map_sp, col=rgb(0,1,0,.05),add=T)# ,lwd=.5,lty=5)
#                     plot(map_sp_2, col=rgb(0,1,0,.05),add=T)
#                     plot(inter_sp, col=rgb(1,0,0,.2),add=T)#,lty=2)
#                     
#                 } else  {
#                     plot(map_sp_2,col=rgb(0,1,0,.05),add=T)
#                     plot(union_comb, col= rgb(0,1,0,.05),add=T,lwd=2)
#                     plot(inter_comb, col=rgb(1,0,0,.2),add=T)#,lwd=2)
#                 }
#             }
#             
#             ##____________ TAB SEQUENCE comparisons  - output  __________### 
#             
#             print(paste(nome[sp],' vs. ', nome[sp_2], ': ',round(D_df[sp,sp_2],3),sep=''))
#             
#             tab_seq = data.frame(rbind(
#                 data.frame(spec = nome[sp], spec_2 = nome [sp_2], spec_ind = sp,
#                            spec_ind_2 = sp_2, D_sim = D_df[sp,sp_2], ratio_IU_cum = ratio_inter_union),
#                 data.frame(spec = nome[sp_2], spec_2 = nome [sp], spec_ind = sp_2, spec_ind_2 = sp, 
#                            D_sim = D_df[sp,sp_2], ratio_IU_cum = ratio_inter_union)))
#             tab_seq = tab_seq %>% mutate_if(is.numeric,funs(round(.,3)))
#             rownames(tab_seq)= NULL
#             if (test_mode) tab_seq
#             
#             if (n_loop == 1) {output[[pattern_out]][['comparison_sequence']] = tab_seq} else {
#                 output[[pattern_out]][['comparison_sequence']] =
#                     data.frame(unique(rbind(output[[pattern_out]][['comparison_sequence']],tab_seq)))    }
#             
#             if (test_mode) output[[pattern_out]][['comparison_sequence']]
#             # mark already analysed (current pair index = lines in df where sp and sp_2 are compared)
#             current_pair_index =  with( output[[pattern_out]][['df']], 
#                                         which( (spec_ind == sp & spec_ind_2 == sp_2) | (spec_ind == sp_2 & spec_ind_2 == sp)))# & output[[pattern_out]][['df']]$spec_ind_2 == sp_2)]
#             # update output
#             output[[pattern_out]][['df']][as.numeric(current_pair_index),'used'] = n_loop
#             
#             if (test_mode) output[[pattern_out]][['df']][as.numeric(current_pair_index),]
#             
#             ##________ NEXT SPP  CHOICE _____________### 
#             # select unused combinations of sp/sp_2
#             d = output[[pattern_out]][['df']] [with(output[[pattern_out]][['df']],
#                                                     which(D_sim > min_sim_vect & used == 0)),] # only unused pairs above D_min are selected / may be redundant min_sim_vect as every D_sim in matrix is above D_sim...
#             # vector with already compared for check 
#             s = output[[pattern_out]][['comparison_sequence']] 
#             D_table_loop = unique(d[order(d[,'D_sim'],decreasing = T),])
#             nxt = 1
#             if (test_mode)d;D_table_loop
#             if (nrow (D_table_loop) < 1) {break; print('next sp pattern: no more combinations above min_sim')}
#             
#             ##____________ NEW SPECIES PAIR TO NEXT INNER LOOP ______##
#             new_sp = D_table_loop$spec_ind[nxt]
#             new_sp_2 = D_table_loop$spec_ind_2[nxt]
#             
#             if (is.na(new_sp) | is.na(new_sp_2)) {break; print('next spp choice: is NA new_sp or new_sp_2')}
#             
#             if (test_mode) paste('new_sp: ',new_sp, 'new_sp_2= ', new_sp_2)
#             
#             ###_______ LOOP to find NEW PAIR in case of repeated new pairs  ________##
#             
#             cond_loop_new_pair = any( new_sp == s$spec_ind & new_sp_2 == s$spec_ind_2 ) | 
#                 any( new_sp_2 == s$spec_ind & new_sp == s$spec_ind_2 )
#             if (test_mode) cat('new pair repeated =', cond_loop_new_pair)
#             while (cond_loop_new_pair){   #if this 'new' pair was already analysed
#                 repeated_pair_index =  with( output[[pattern_out]][['df']], 
#                                              which( (spec_ind == new_sp & spec_ind_2 == new_sp_2) |
#                                                         (spec_ind == new_sp_2 & spec_ind_2 == new_sp)))# & output[[pattern_out]][['df']]$spec_ind_2 == sp_2)]
#                 output[[pattern_out]][['df']][as.numeric(repeated_pair_index),'used'] = n_loop
#                 nxt = nxt + 1
#                 if (nxt > length (D_table_loop)) {loop = 'break'; cat(
#                     'loop for repeated new pairs: no more combinations above min_sim  \n')
#                 break
#                 }
#                 
#                 # Trying with a new species pair
#                 new_sp = D_table_loop$spec_ind[nxt]
#                 new_sp_2 = D_table_loop$spec_ind_2[nxt]
#                 # paste(nome[new_sp],new_sp,' \ ', nome[new_sp_2], new_sp_2, 'D=', D_df[new_sp,new_sp_2])
#                 #
#             } # < - WHILE  new pair END
#             
#             ##_______new_sp / new_sp_2 to sp / sp_2  assignment _______##
#             if (loop == 'break') {  break ; cat('loop for repeated new pairs \n new_sp to sp: no more combinations above min_sim') 
#                 # assign new species to sp and sp_2
#             }  else  { sp = new_sp; sp_2 = new_sp_2; D_loop = D_df[sp,sp_2]    }
#             
#             # announce new comparison inside pattern
#             cat ('new comparison: ')
#             cat(nome[new_sp],new_sp,' X ', nome[new_sp_2], new_sp_2, '; D =', D_df[new_sp,new_sp_2],'\n')
#             cat('pattern_out = ',pattern_out,'/ ratio_inter_union = ',ratio_inter_union,'/ D_loop = ', D_loop )    
#             # output;sp;sp_2
#             paste(nome[new_sp],new_sp,' X ', nome[new_sp_2], new_sp_2, '; D =', round(D_df[new_sp,new_sp_2],4), '/ loop =', n_loop)
#             
#             ##______ second BREAK AsSEssMENT  ________________##
#             # just to check loop manualy # check loop conditions
#             if ((D_df[sp,sp_2] < min_sim) | is.na(D_df[sp,sp_2]) ) { break; cat('second BREAK AsSEsMENT: D_loop sp/sp_2 below min_sim') } # meaning there is no sp_2 with signficant D_sim to sp_pattern; folows the next species in for loop
#             if ((loop == 'break') | (D_loop < min_sim) | (ratio_inter_union < ratio_IU)) {break; cat('second BREAK AsSEsMENT: (loop == break) | (D_loop < min_sim) | (ratio_inter_union < ratio_IU) ') }
#             if (sp != sp_pattern) { 
#                 # HERE check only after loop and if sp != sp_pattern 
#                 if ( (D_df[sp_pattern,sp] < min_sim_to_pattern) |
#                      (D_df[sp_pattern,sp_2] < min_sim_to_pattern) |   # Dist to original species must be above min_sim_to_pattern
#                      (is.na(D_df[sp_pattern,sp])) |
#                      (is.na(D_df[sp_pattern,sp_2]))) { 
#                     if (spatial_check) { check_spat_results = check_spatial()
#                     paste('check_spat_results = ',check_spat_results[['status']])
#                     if (check_spat_results[['status']]) {print('second BREAK AsSEsMENT: spatial test OK ') } else {  break; print('D_sim with bad value but spatial test D<min_sim ' ) }     #if there is problem with dist to pattern, check must be done spatially...
#                     } #spatial_check ends
#                 }
#             } #if sp!= sp_pattern ends  
#             #__________________________________________________________________________\
#         }  # end internal loop repeat   ####   inside pattern                   )
#         #_________________________________________________________________________/    
#         ################   EXTERNAL LOOP - new pattern assignement
#         if (test_mode) cat ('mark_1 - end internal loop')
#         # species pool of the whole pattern
#         df = output[[pattern_out]][['comparison_sequence']]
#         spp_pool_pattern = unique(unique(df$spec),unique(df$spec_2))
#         
#         ################   LAST SPATIAL CHECK    
#         if (spatial_check)  {  # need an apropriate place because only the last round will be computed (??june2018??) - Neverthless it is in inner loop as sp sp-2 changes..
#             output[[pattern_out]][['summary']] = data.frame(D_init, ratio_inter_union, D_loop, check_spat_results[['D_chk']], check_spat_results[['D_chk_sp_pattern']],
#                                                             check_spat_results[['D_chk_sp_2_pattern']], as.character(check_spat_results[['status']]))
#             names(output[[pattern_out]][['summary']]) = c('D_init', 'ratio_inter_union', 'last D_loop',
#                                                           'last_D_chk', 'last_D_chk_sp_patt','last D_chk_sp_2_patt','last check_spat_results')
#             #round in a tidy format
#             if (test_mode) output[[pattern_out]]['summary']
#             if (test_mode) cat('mark_2 - output writen')
#         }
#         ###_____   SAVING union / intersection _________________##
#         if (save_spat_obj){
#             # if (!is.na(inter_comb))  output[[pattern_out]][['intersection']]=inter_comb
#             # if (!is.na(union_comb))  output[[pattern_out]][['union']] = union_comb   
#             dir_shp = "C:/SIG2018/CO_OCCUR_PATT/patterns_shp"
#             layer_shp_inter = paste( pattern_out,'_inter', sep = '')
#             layer_shp_union = paste( pattern_out,'_union', sep = '')
#             st_write(inter_comb, dsn = dir_shp, layer = layer_shp_inter, driver = "ESRI Shapefile", update = TRUE)
#             st_write(union_comb, dsn = dir_shp, layer = layer_shp_union, driver = "ESRI Shapefile", update = TRUE)
#             if (test_mode) cat("mark_3 saved shp")
#         }
#         
#         cat(pattern_out," Ends",'\n') 
#         
#         # write csv with table of species inside pattern
#         output_name = paste (dir_output, pattern_out, '_output.csv', sep = "")
#         write.csv2(output[[pattern_out]][['df']],output_name)
#         if(pdf_plot) { plot_file_name = paste(pattern_out,'%d.pdf',sep = "")
#         dev.copy(pdf, plot_file_name)#, width=15, height=15*.6)
#         dev.off()
#         n_dev_copy = 0
#         }
#         if (test_mode) cat('mark_4 - end external (sp) loop')
#         ###_______________________________________________________________________________________\
#         ###________ END EXTERNAL loop for  _________________________________________________________##                                      
#     }
#     ###_______________________________________________________________________________________|
#     options(warn = oldw)
#     return(output)
#     } 
# # end function co_occ and sugestions

#___________________________________________________________________________________________________________________________#
# plot_spp_focus 27NOV2018
'plot sp focus in range sim max-min, choose between internal_to, overlap and contains in position parameter'
'Saved objects =  df_spp species list used to print'
'color_scheme = 1, 2 or 3 are color gradients; 4 is all blue and other options are direct, e.g. color_scheme = rgb(.5,.3,.7,.2)'
'methods = "sp_focus" or "spp_df" in which the data frame must have SCINAME, C_thresh, status_sp_focus & area column names'
plot_spp_focus = function ( method = 'sp_focus', # option method = 'spp_df'
                            sp_focus,  spp_df = spp_df, max_sim = 1, min_sim = .1, similarity_df = sim_df,
                            area_df = area_species,  position = FALSE, # plots everything and # position = c('contains', 'overlap') selects internal spp
                            map_species = mapAE_SA, backgroud_map = sa, new_window = FALSE, new_plot = TRUE, color_scheme = 1,...)  { # color_scheme = rgb(0,0,1,.005)
    # area_species = sim_df_lowlands %>% select(SCINAME = species1,area = area_sp1) %>% unique (can use to calculate areas)
    # print('plot_spp_focus_V1.1')#color_scheme
    print('plot_spp_focus_V1.2')#color_scheme == 4 (all blue), method = spp_df, etc 27/NOV/2018
    print(paste('method =',method))
    if (!exists("area_species")){area_df = area_species_SA}
    if (method != 'sp_focus' & method != 'spp_df'){ print('method must be specified: "sp_focus" or "spp_df"');break }
    if (method == 'sp_focus'){
        sp_df = sim_sp_df(sp_focus = sp_focus, max_sim = max_sim ,min_sim = min_sim)
        
        # int_ext parameters from plot_spp_focus
        sp_df = sp_int_ext(sp_focus = sp_focus, tab_spp = sp_df) # sp_int_ext to verify a df to position internal, external or overlay
        
        print(sp_focus)
        print(paste('N-spp: ',length(sp_df$SCINAME),sep=''))
        sp_df %>% group_by(status_sp_focus) %>% summarise(n=n()) %>% print
        #sort by area
        if(!is.character(area_df[,'SCINAME'])){area_df = area_df %>% mutate_if(is.factor,as.character)}
        if(!is.numeric(area_df[,'area'])){area_df[,'area'] = area_df %>% select(area) %>% mutate_if(is.character,as.double)}
        sp_df = sp_df %>% left_join(area_df, by=c('SCINAME'='SCINAME')) %>% arrange(desc(area))
        # filter species overlap or contains
        if(position[1] != FALSE) {sp_df = sp_df %>% filter(status_sp_focus %in% position) %>% as.data.frame() }
    }
    if (method == 'spp_df'){
        if (!exists('sp_focus')) {sp_focus = 'species focus'}
        if (all(c('SCINAME','C_thresh','status_sp_focus','area') %in% names(spp_df))){ sp_df = spp_df;
        sp_df = sp_df %>% arrange(desc(area))
        } else { print('columns spp_df must be "SCINAME","C_thresh","status_sp_focus" & "area"')}
    }
    #sort by area
    
    # plot 
    if(new_window) {x11()}
    bb_max = bb_max_spp(sp_df, map_base = map_species)
    if(new_plot) {
        title = paste( paste(str_sub(sp_focus,1,1),'.', word(sp_focus,2,2), sep = ''),' N=',nrow(sp_df),' ',round(max_sim,2),'-',round(min_sim,2), sep='')
        plot(st_geometry(backgroud_map), main = title,
             col = rgb(.7,.6,.6,.5), xlim = bb_max[c(1,3)], ylim = bb_max[c(2,4)])}
    if (color_scheme > 3){ color_scheme = rgb(0,0,1,.005)}
    
    for (i in 1:length(sp_df$SCINAME)){
        graph_factor = ifelse(i==1,.1,((i-1)/(length(sp_df$SCINAME)-1)))
        map = map_species %>% filter(SCINAME == sp_df$SCINAME[i]) %>% group_by(SCINAME) %>% summarise() %>% st_geometry
        if (color_scheme == 1 | color_scheme == 2 | color_scheme == 3){
            if(color_scheme == 1){color_scheme = rgb(1-graph_factor/1.4,1-graph_factor/2,graph_factor/1.7, graph_factor/2 )
            } else { if(color_scheme == 2) {color_scheme = rgb(1-graph_factor/1.2,graph_factor/1.5,1-graph_factor/1.1, graph_factor/1.8 )
            } else { color_scheme = rgb(graph_factor/1.8,graph_factor/1.3,1-graph_factor/1.7, graph_factor/2 ) }}
        }
        
        plot(map, col = color_scheme, add = T)
    }
    plot(map_species %>% filter(SCINAME == sp_focus) %>%
             group_by(SCINAME) %>% summarise() %>% st_geometry,lwd=3,add=T)
    sp_df
}

if(examples){plot_spp_focus (sp_focus = 'Capito auratus', min_sim = .7, new_window = F, position = FALSE)
    plot_spp_focus (sp_focus = 'Capito auratus', max_sim=.6, min_sim = .5, new_window = F,
                    position = c('contains','overlap',NA))}
# results of the function
# 1 contains            3
# 2 internal_to         3
# 3 NA                 14
# [1] "bb_max_spp V1.1"
# sp_focus status_sp_focus                    SCINAME C_thresh IU_sp_sp area_sp_sp_ratio         area
# 1  Capito auratus     internal_to          Penelope jacquacu   0.6150    0.617            0.623 5.007343e+12
# 2  Capito auratus     internal_to      Ochthornis littoralis   0.6474    0.653            0.684 4.562556e+12
# 3  Capito auratus            <NA>       Automolus infuscatus   0.6132    0.630            0.721 4.323531e+12
# ...
# 18 Capito auratus        contains Epinecrophylla haematonota   0.6758    0.680            1.416 2.201917e+12


# area_df = area_species_SA, position = c('contains', 'overlap'), position = FALSE plots everything
# map_species = mapAE_SA, backgroud_map = sa, new_window = FALSE, add_plot = FALSE)
#___________________________________________________________________________________________________________________________#
# plot_spp_sim_maps -> spp on a vector function 
{'PLOT SPECIES in a df (SCINAME, dist-sim) and calculates intersection and union as result'
    'map_species1 = map with all species in a SCINAME col & background_map1 = neotropical geometry'
    'if plot_level_inter_union the function coherence will plot just the level inter union, not accummulated across levels'
    'RESULTS in a list with UNION AND INTERSECTION'}
# tab_spp = plot_spp_focus (sp_focus = 'Capito auratus', min_sim = .6, new_window = F, position = FALSE, sim_df = sim_df_lowlands, area_df = area_species) %>% select(SCINAME,C_thresh)}
plot_spp_sim_maps = function (tab_spp = tab_spp,
                              map_species = mapAE_SA,
                              n_levels = FALSE,
                              background_plot = FALSE, # default option to use with coherence_to_sp function!
                              background_map = sa, # if background_map = FALSE no background is ploted
                              title = FALSE, # if title != FALSE needs to be given as character
                              plot_level_inter_union = FALSE,... ) {
    print('this is plot_spp_sim_maps V1')
    tab_spp = tab_spp %>% select('SCINAME'= 1, 'C_thresh'=2) %>% arrange(desc(C_thresh)) %>%
        mutate_if(is.factor,as.character)%>% as.data.frame # not working as tibble!!!
    ## PLOT background
    if (background_plot) {
        if(title != FALSE) {plot (background_map$geometry, main = title, cex.main = 1.3)
        } else {plot (background_map$geometry)}
    }
    
    if (n_levels != FALSE){ group_lvl = ntile(tab_spp$C_thresh,n_levels) }
    
    ## PLOT all species in the level
    for(i in 1 : nrow (tab_spp)) {
        m = map_species %>% filter(SCINAME == tab_spp[i,'SCINAME']) %>% st_geometry
        plot(m, add=T, col=rgb(0,1,0,0.2))
        # intersection # union
        if (i == 1){ inter = m; union = m; next }
        inter_backup = inter
        inter = m %>% st_intersection(inter)
        union = m %>% st_union (union)
        if (st_crs(inter) != st_crs(map_species)) {
            st_crs(inter) = st_crs(map_species);
            st_crs(union) = st_crs(map_species)
        }
        if (length (st_area(inter)) == 0) { inter = inter_backup # keeps the old interesection
        print (paste ('BREAK not all species have commom area -->> ', tab_spp[i,'SCINAME']));
        break }
    }
    # plots outside loop
    print('fim_loop')
    if(plot_level_inter_union){
        # plot (union, add = T, col = rgb (0,0,1,.15), lwd = 2, lty = 1)
        # plot (inter, add = T, lwd = 2, col = rgb (.7,.5,0.2,.8))}#(0,.8,.8,.8)) # the last non-zero interception
        plot(st_geometry(union), add = T, col = rgb(0,1,1,.01), lwd = 2, lty = 1)
        plot(st_geometry(inter), add = T, col = rgb(1,.5,.3,.7), lwd = 2)
    }
    list('union' = union, 'intersection' = inter)
}


#___________________________________________________________________________________________________________________________#
# plot_one_by_one - graphs comparing sp focus with each species in a table
'spp is a data.frame with SCINAME, status_sp_focus and C_thresh columns'
'adjust the row_col parameter to the number of total comparisons e.g. 10 maps ~ row_col = c(2,5)'
'V9.5 map to map_species'

plot_one_by_one = function (spp = g.ruf_ext, sp_focus = 'Gymnopithys rufigula', map_species = mapAE_SA,
                            row_col = c(2,5), new_window = TRUE, col1 = rgb(0,1,0,.2), col2 = rgb(0,0,1,.2), ...){
    spp %>% nrow %>% print # check number of species and set mfrow
    spp = spp %>% arrange(status_sp_focus,desc(C_thresh))
    if(new_window) {x11(); par(mfrow = row_col)}
    sp_map = map_species %>% filter(SCINAME == sp_focus)
    for(i in 1:nrow(spp)) {
        m = map_species %>% filter(SCINAME == spp$SCINAME[i])
        plot(st_geometry(m),col = col1, main = paste(abbrev(spp$SCINAME[i]), spp$C_thresh[i], spp$status_sp_focus[i]))
        plot(st_geometry(sp_map), col = col2, add=T)
    }}


#________________________________________________________________________________________________________________________________________________
'calculates inter union ratio after incorporating each species sorted by sim'
'results in list with [[df]] [[intersection]] [[union]]; df has the ratio inter_union accumulated and areas inter union for each combination'

inter_union_accum = function (sp_focus, spp_df, map_species = mapAE_SA,
                              plot_map = FALSE, plot_add = TRUE,
                              ...) {
    print('inter_union_accum V1') # 29NOV2018
    df = spp_df %>% arrange(desc(C_thresh))
    sp_focus_map = map_species %>% filter(SCINAME == sp_focus)
    #if(plot_map){ plot(st_geometry(sp_focus_map),col = rgb(0,0,1,.7), main = paste(abbrev(sp_focus),'N=',nrow(df))) }
    intersection = union = sp_focus_map
    df_res = data_frame('IU_accum' = double(), 'area_inter' = double(),'area_union' = double())
    for (i in 1:nrow(df)){
        sp_map = map_species %>% filter(SCINAME == df[i,'SCINAME'])
        intersection = st_intersection(sp_map, intersection)
        union = st_union(sp_map, union)
        intersection_area = st_area(intersection) %>% as.double
        union_area = st_area(union) %>%  as.double
        df_res[i,] = c(round(as.double(intersection_area/union_area),3), as.double(intersection_area), as.double(union_area))
        if(plot_map){  plot(union, col = rgb(0,1,0,1/nrow(df)), add = plot_add)
            plot(intersection, col = rgb(1,0,0,1/nrow(df)), add = plot_add)   
        }
        if(plot_map){ plot(st_geometry(union), lwd = 2.5, add = plot_add)
            plot(st_geometry(sp_focus_map), lwd = 3, add = plot_add)
            plot(st_geometry(intersection), lwd = 2.5, add = plot_add)}
    }
    result = list('df' = cbind(df,df_res), 'intersection' = intersection, 'union' = union)
}

# sugestions: ?
if(examples){
    test = inter_union_acumm(sp_focus = 'Gymnopithys rufigula', spp_df = g.ruf_int[1:20,], map_species = mapAE_SA, plot_map = TRUE, plot_add = TRUE)
    test[['df']] }
#________________________________________________________________________________________________________________________________________________


#___________________________________________________________________________________________________________________________#
'INTER_PATT - plot & calculate & save shape of intersection and union of a group of species'
# 'used in sp_levels_sim for species in the "level"'
'relies on "bb_max_spp"'

inter_patt = function (sp_df = species_df, sort_by_sim = FALSE, title = 'title', map_species = mapAE_SA, plot_background = TRUE, new_window = FALSE, background_map = sa,... ){
    # print('inter_patt V1.1 29NOV2018') # option NO plot; mapSA internal was exchanged to call object 'map_species' using mapAE_SA in default
    print('inter_patt V1.2 17MAI2019') # map_species in bb_max_spp call
    break_function = FALSE
    
    bb_max_spp_lvl = bb_max_spp (tab_spp = sp_df, column_spp = 'SCINAME', map_base = map_species)
    
    # check for NAs caused by spelling
    if(bb_max_spp_lvl %>% is.na %>% any) {
        sp_df[,column_spp]  = sp_df[,column_spp] %>% gsub('_', ' ',.)
        bb_max_spp_lvl = bb_max_spp (tab_spp = sp_df, column_spp = 'SCINAME', map_base = map_species)
    }
    # check again changing space to underline
    if(bb_max_spp_lvl %>% is.na %>% any) {
        sp_df[,column_spp]  = sp_df[,column_spp] %>% gsub(' ', '_',.)
        bb_max_spp_lvl = bb_max_spp (tab_spp = sp_df, column_spp = 'SCINAME', map_base = map_species)
    }
    
    if (sort_by_sim == TRUE & 'C_thresh' %in% names(sp_df)) { sp_df = sp_df %>% arrange(desc(C_thresh)) }
    map_1 = map_species %>% filter(SCINAME == sp_df[1,'SCINAME']) %>% st_geometry
    
    if(plot_background){
        if(new_window) { x11() }                #; par(mar=c(1,3,1,1)) }
        plot(background_map$geometry, xlim = bb_max_spp_lvl[c(1,3)], ylim = bb_max_spp_lvl[c(2,4)], col=rgb(0,.3,0,.4))
        title(main = title, cex.main = 2)
    }
    plot(map_1, add = T, col = rgb(1,0,0,.1))
    
    for(i in 2:nrow(sp_df)){
        m = map_species %>% filter(SCINAME == sp_df[i,'SCINAME']) %>% st_geometry
        plot(m, add = T, col = rgb(0,0,1,0.1))
        if (i ==2) {union = m %>% st_union(map_1)} else { union = m %>% st_union(union) %>% st_geometry() }
        if (i ==2) {inter = m %>% st_intersection(map_1) %>% st_geometry} else { inter = m %>% st_intersection(inter) %>% st_geometry }
        
        if (length(st_area(inter)) == 0) {
            print('not all species have commom area')  # not sure if this break is working
            print(paste(sp_df[i,'SCINAME'], 'in', sp_df[,'SCINAME']))
            break_function = T
            break }
        
        
        # plot(inter, add = T,col = rgb(0,1,0,( i/nrow( sp_df ))))
    }
    if(break_function) {print('No intersection till the end')
        plot(background_map$geometry, xlim = bb_max_spp_lvl[c(1,3)], ylim = bb_max_spp_lvl[c(2,4)], col=rgb(0,.3,0,.4))
        break }
    
    plot(inter, add = T, col = rgb(0,1,1,1),lwd = 1.4) # red the last non-zero interception
    plot(union,add = T, col = rgb(.4,1,1,.02), lwd = 2)
    
    maps = list()
    maps[['intersect']] = inter  # the returned object is the intersection sf shape
    maps[['union']] = union
    
    maps      # list with intersection and union 
}
#___________________________________________________________________________________________________________________________#

#___________________________________________________________________________________________________________________________#
# ex-tab_spp wraper to make fast tab_spp tables [,c(SCINAME, C_thresh)]
sim_sp_df = function (sp_focus, max_sim = 1, min_sim = .7, similarity_df = sim_df) {
    similarity_df %>% filter(species1 == sp_focus) %>% filter (species2!= sp_focus) %>% 
        select('SCINAME' = species2, C_thresh) %>% filter(C_thresh >= min_sim & C_thresh <= max_sim) %>%
        mutate_if(is.factor, as.character) %>% arrange(desc(C_thresh)) %>% as.data.frame()
}
#__


# Smaller tools
# BB MAX SPECIES 
'for vector of species or data frame'
bb_max_spp = function (tab_spp = tab_spp, column_spp = 'SCINAME', map_base = mapAE_SA) { 
    print('bb_max_spp V1.1') # tab_spp is a df with species in a column SCINAME; map_base object as a parameter
    box_df = data_frame ()
    if (is.vector(tab_spp)){ tab_spp = data.frame( 'SCINAME' = tab_spp, stringsAsFactors=F)} # check column_spp, if matches the name, ok
    # tab_spp = tab_spp %>% mutate_if(is.factor,as.character)
    for (i in 1:nrow(tab_spp)){
        # map & data frame with all boxes
        map1 = map_base %>% filter(SCINAME == tab_spp[i,column_spp]) %>% st_geometry
        box_df = rbind(box_df,data_frame('xmin'= st_bbox(map1)[1],'ymin'= st_bbox(map1)[2],
                                         'xmax'= st_bbox(map1)[3],'ymax'=st_bbox(map1)[4] ))
    }
    #check spelling '_' vs. ' '
    if ( box_df %>% is.na %>% all) { 
        box_df = data_frame ()
        tab_spp[,column_spp] = tab_spp[,column_spp] %>% gsub('_', ' ',.)
        for (i in 1:nrow(tab_spp)){
            # map & data frame with all boxes
            map1 = map_base %>% filter(SCINAME == tab_spp[i,column_spp]) %>% st_geometry
            box_df = rbind(box_df,data_frame('xmin'= st_bbox(map1)[1],'ymin'= st_bbox(map1)[2],
                                             'xmax'= st_bbox(map1)[3],'ymax'=st_bbox(map1)[4] ))
        }
    }
    
    xy = c('xmin'=min(box_df$xmin), 'ymin'=min(box_df$ymin), 'xmax'=max(box_df$xmax), 'ymax'=max(box_df$ymax))
    xy
}
if(examples) { 
    bb_max_test = bb_max_spp(tab_spp, map_base = mapAE_SA)}
#___________________________________________________________________________________________________________________________

# check lots of species in a pattern against a list of species patterns
check_spp_patt = function(sp_patt, spp_chk, thr_df = thr, only_name_number = T) {
    for (sp in spp_chk){
        sp_patt_list = thr_df %>% filter(species1 == sp_patt) %>% group_by(species1,species2) %>%
            summarise(n = n(), sim = mean(C_thresh), min = min(min), max = max(max)) %>%
            group_by(species2, n, sim, min, max) %>% summarise()
        df = thr_df %>% filter(species1 == sp ) %>% group_by(species1,species2) %>% 
            summarise(n = n(), sim = mean(C_thresh), min = min(min), max = max(max)) %>% 
            group_by(species2, n, sim, min,max) %>% summarise()    
        n = inner_join(df, sp_patt_list, by = c('species2'))
        if(only_name_number) { print(paste(sp_patt,'    x     ',sp,'   :',  nrow(n), ' spp  ', ifelse(nrow(n)==0,'...ZERO!!!!','')))} else {print(n)}
    }
}
#_____________________________________________________________________________________________________________________________


# wrap to make fast list of sim_df 
list_df = function (sp_focus, max_sim, min_sim, df_sim = sim_df){
    print('list_df V2')# option to sp_focus
    if(exists('sp_focus')){
        df = df_sim %>% filter(species1 == sp_focus) %>% filter(species1 != species2) %>%
            filter (C_thresh <= max_sim & C_thresh >= min_sim) %>% arrange(desc(C_thresh))
    } else {
        df = df_sim %>% filter(species1 != species2) %>%
            filter (C_thresh <= max_sim & C_thresh >= min_sim) %>% arrange(desc(C_thresh))}
    df
}

# _________________________________________________________________________
abbrev = function(x){paste(str_sub(x,1,1),'.', word(x,2,2), sep = '')}

#_____________________________________________________________________________
ggmap_bbox <- function(map) {
    if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
    # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
    # and set the names to what sf::st_bbox expects:
    map_bbox <- setNames(unlist(attr(map, "bb")), 
                         c("ymin", "xmin", "ymax", "xmax"))
    
    # Coonvert the bbox to an sf polygon, transform it to 3857, 
    # and convert back to a bbox (convoluted, but it works)
    bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))
    
    # Overwrite the bbox of the ggmap object with the transformed coordinates 
    attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
    attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
    attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
    attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
    map
}
if(examples){
    # Use the function:
    # map <- ggmap_bbox(map)
    # ggmap(map) + 
    #     coord_sf(crs = st_crs(3857)) + # force the ggplot2 map to be in 3857
    #     geom_sf(data = nc_3857, aes(fill = AREA), inherit.aes = FALSE)
    #     
    map_g.rufigula = mapAE_SA %>% filter(SCINAME =='Gymnopithys rufigula') %>% group_by(SCINAME) %>% summarise()
    g.rufi_3857 = st_transform(map_g.rufigula, crs = 3857)
    # using plot and bgMap
    plot(g.rufi_3857, bgMap = sa_ggmap_z4,col=rgb(0,0,1,.1),lwd=2)
    # Using the function to modify the bbox of ggmap obj to 3857 
    sa_ggmap_z4_ggbb <- ggmap_bbox(sa_ggmap_z4)
    ggmap(sa_ggmap_z4_ggbb) + 
        coord_sf(crs = st_crs(3857)) + # force the ggplot2 map to be in 3857
        geom_sf(data = st_transform(map_g.rufigula,crs = 3857), aes(), alpha = .2, col='blue',inherit.aes = FALSE)
} 
# ____________________________________________________________________________________________________________________________________________

# Check overlap in species from patterns (sp_patt is a vector with one spp; spp_chk vector with many spp)
# uses thr data frame
check_spp_patt = function(sp_patt, spp_chk, thr_df = thr, only_name_number = T) {
    for (sp in spp_chk){
        sp_patt_list = thr_df %>% filter(species1 == sp_patt) %>% group_by(species1,species2) %>%
            summarise(n = n(), sim = mean(C_thresh), min = min(min), max = max(max)) %>%
            group_by(species2, n, sim, min, max) %>% summarise()
        df = thr_df %>% filter(species1 == sp ) %>% group_by(species1,species2) %>% 
            summarise(n = n(), sim = mean(C_thresh), min = min(min), max = max(max)) %>% 
            group_by(species2, n, sim, min,max) %>% summarise()    
        n = inner_join(df, sp_patt_list, by = c('species2'))
        if(only_name_number) { print(paste(sp_patt,'    x     ',sp,'   :',  nrow(n), ' spp  ', ifelse(nrow(n)==0,'...ZERO!!!!','')))} else {print(n)}
    }
}

## THIS VERSION WORKS AFTER source_V9.5
# Check overlap in species from patterns (sp_patt is a vector with one spp; spp_chk vector with many spp)
# uses thr1 data frame for varzea data frame structure (species1,species2,Csim,Cthres_max,Cthres_min, depth_min, depth_max)
check_spp_patt_2 = function(sp_patt, spp_chk, thr_df = thr, only_name_number = T) {
    results_df = data_frame()
    for (sp in spp_chk){
        sp_patt_list = thr_df %>% filter(species1 == sp_patt) %>% group_by(species1,species2) %>%
            summarise(n = n(), Cs_mean = mean(Csim), Ct_min = min(Cthres_min), Ct_max = max(Cthres_max)) %>%
            group_by(species2, Cs_mean, Ct_min, Ct_max) %>% summarise(n())
        
        df = thr_df %>% filter(species1 == sp ) %>% group_by(species1,species2) %>% 
            summarise(n = n(), Cs_mean = mean(Csim), Ct_min = min(Cthres_min), Ct_max = max(Cthres_max)) %>% 
            group_by(species2, Cs_mean, Ct_min, Ct_max) %>% summarise()
        
        n = inner_join(df, sp_patt_list, by = c('species2')) %>% select(common_spp = species2, everything())
        N_spp_patt = thr_df %>% filter(species1 == sp_patt) %>% group_by(species1,species2) %>% nrow()
        N_spp_chk = thr_df %>% filter(species1 == sp) %>% group_by(species1,species2) %>% nrow()
        if(only_name_number) { 
            # print(paste(sp_patt, '(',N_spp_patt,'spp)', '    ',  nrow(n), ' spp', '     ', sp, '(',N_spp_chk,'spp)', ifelse(nrow(n)==0,'...ZERO!!!!','')))
            results_df = rbind(results_df, data_frame(species1 = sp_patt, nspp1 = N_spp_patt, common_spp = nrow(n), nspp2 = N_spp_chk, species2 = sp))%>%
                arrange(desc(common_spp), desc(nspp2))
        } else {print(sp); print(n)}
    }
    results_df
}

# print('This is source_clean_V9.5')
print ('This is congruence_V1.0.R')
