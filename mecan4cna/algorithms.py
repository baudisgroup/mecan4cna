
# coding: utf-8

# In[1]:
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from mecan4cna.detect_peaks import detect_peaks
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd



class mecan:
    def __init__(self, bins_per_interval=20, intervals=4, interval_step=2,  minprobes=1000, probe_thresh=3, mpd_coef=0.1,
     outpath=None, plot=False, verbose=False, ranking_method='ass', thresh_min=0.3, thresh_max=1.3, minlevel=9,
     neglect_level_cap_high = 0.8, neglect_level_cap_low = 0.2, neglect_thresh = 0.3 ):

        self.bins_per_interval = bins_per_interval
        self.intervals = intervals
        self.interval_step = interval_step
        self.minprobes = minprobes
        self.probe_thresh  = probe_thresh
        self.mpd_coef = mpd_coef
        self.outpath = outpath
        self.plot = plot
        self.verbose = verbose
        self.ranking_method = ranking_method
        self.thresh_min=thresh_min
        self.thresh_max=thresh_max
        self.minlevel = minlevel
        self.neglect_level_cap_high = neglect_level_cap_high
        self.neglect_level_cap_low = neglect_level_cap_low
        self.neglect_thresh = neglect_thresh
        self.total_probes = 0





    # Functions



    # Find all signal peaks of a segment file,
    # 
    # LATEST VERSION
    #
    #### Params  ####
    # segment: the input segments, a list of dictionary. keys: chrom, start, end, probes, value
    # bins_per_interval: the number of bins between two level intervals, default=20
    # plot: wheather to generate histogram, default=False
    # minprobes: the minimumn number of probes a bin must have in order to be counted. default=500
    #
    #### Returns ####
    # A dataframe of peaks, columns: value(totoal probe counts), bin(CN level)
    #
    #### Misc ####
    # Hard coded varaiables
    # intervals = 4, 0 - 5 copies
    # probe_thresh = 3, minimumn number of probes of a segment
    def computePeaks(self, segments, bins_per_interval=None, plot=False):

        if bins_per_interval is None:
            bins_per_interval = self.bins_per_interval

        # init
        bin_size = 1/bins_per_interval
        segbin = []
        for i in range(bins_per_interval * self.intervals):
            segbin.append(0)
        
        # loop through all segments
        for seg in segments:
            
            # convert log ratio back to copy number value
            value = round(2**seg['value'] * 2, 4)
            probes = int(seg['probes'])

            # drop segments with too few probes
            if probes >self.probe_thresh:
                # take the log10 of probes as the weight value to put into bins
    #             weighted = np.log(probes) 
                weighted = probes

                # ingore extrem values
                if value < (self.intervals+bin_size/2) and value > bin_size/2:
                    # find in which bin to put the value in
                    # bin index is determined by the VALUE of the segment
                    # the adding up value is the WEIGHTED 
                    bin_index = math.floor( (value-bin_size/2) * bins_per_interval)
    #                 bin_index = math.floor( (value-0.1) * bins_per_interval)
                    # adding up
                    segbin[bin_index] += weighted   
                    
        # filter by number of probes
        for i in range(len(segbin)):
            if segbin[i] < self.minprobes:
                segbin[i] = 0

        #####  peak calculation ####
        
        # bin scales
        bin_scales = np.arange(bin_size, self.intervals+bin_size, bin_size)
        
        # compute mpd as 20% of bins_per_interval
        mpd = round(bins_per_interval * self.mpd_coef)
        peaks_index = detect_peaks(segbin, mpd=mpd)
        peaks_value = [segbin[i] for i in peaks_index]
        peaks_bin = [bin_scales[i] for i in peaks_index]
        
        

        # plot
        if plot or self.outpath:
            plt.figure(figsize=(30,5))
            plt.xticks(rotation=70)
            plt.bar(bin_scales, segbin, width=bin_size/2, tick_label= np.round(bin_scales, 2))
            plt.plot(bin_scales, segbin, 'orange')
            plt.plot(peaks_bin, peaks_value,'rx')
            plt.xlabel('Virtual Copy Number Levels')
            plt.ylabel('Number of Probes')

            if self.outpath:
                try:
                    plt.savefig(os.path.join(self.outpath, 'histogram.png'), bbox_inches='tight')
                except Exception as e:
                    print(e)
            if plot:
                plt.show()

                 


        df_peaks = pd.DataFrame({'value': peaks_value, 'bin': peaks_bin})
        df_peaks['bin'] = round(df_peaks['bin'],2)   
        return df_peaks




    def levelScore(self, peaks, base, thresh ):
        levels = peaks.loc[:,['bin', 'value']]
        levels['R_level'] = round((levels['bin'] - base) / abs(base-thresh)).astype(int)
        dup_sum = sum(levels[levels.R_level >0].R_level.drop_duplicates())
        del_sum = sum(levels[levels.R_level <0].R_level.drop_duplicates())
        # return [sum(levels.R_level), dup_sum, del_sum]
        return [dup_sum+del_sum, dup_sum, del_sum]


    # Choose one peak as the baseline, another peak as standard distance (one copy difference)
    # Then, evaluate the distance of other peaks to the closest integral multiple of this distance to the base.
    #
    # Mainly used in functions computeModels()
    #### Params ####
    # df: input dataframe, expecting output of function computePeaks()
    # bindex: index of the baseline
    # lindex: index of standard level
    #
    #### Returns ####
    # the sum of the scaled offsets from integral multiples 
    # Note: the scaling to remove the bias from number of probes. 
    # it reflects the weight of bins with high probes.
    def modelScore(self, df, bindex, lindex, showtable=False):
        # compute distance to the baseline, for each peak
        df['dist'] = abs(df['bin'] - df.loc[bindex, 'bin'])
        
        # divide by the standard distance
        df['cn'] = df['dist'] / df.loc[lindex, 'dist']
        
        # round to an integer as the integral multiple
        df['round'] = round(df['cn'])
        
        # compute the difference between the real and rounded multiple.
        df['off'] = round(abs(df['cn'] - df['round']), 2)
        
        # the max number of probes except the base and starndard level
    #     max_value = max(df['value'].drop(df['value'].index[[bindex,lindex]]))
        max_value = max(df['value'])
    #     # weight = df['value'] / max_value
    #     df['scaled_off'] = df['off'] * df['value'] / max_value
        df['weight'] = df['value']/max_value
        df['scaled_off'] = df['off'] * df['weight']
    #     off_sum = sum(df['scaled_off'])
    #     off_sum = off_sum * (df.loc)
        
        df_neg = df[(df.cn < self.neglect_level_cap_high) & (df.cn > self.neglect_level_cap_low)]


        if showtable:
            print(df)
        return sum(df['scaled_off']), sum(df_neg['value'])
    #     return sum(df['off'])

    # Compute model scores of all the possible combinations of baseline and standard levels
    # Only half of the matrix need to be computed, the other half is the same, and diagnal doesn't make sense.
    #
    #### Params ####
    # df_peaks: dataframe of detected peaks, output of function computePeaks()
    # peaktable: wheather to print the peak table, default= False
    #### Returns ####
    # A dataframe with all the combinations and scores
    def computeModels(self, df_peaks, showtable=False, lcap=11, slcap=7):

        # create an empty table
        peak_table = [[None] * df_peaks.shape[0] for i in range(df_peaks.shape[0])]
        neg_table = [[None] * df_peaks.shape[0] for i in range(df_peaks.shape[0])]

        # compute all possible combinations and fill the table
        for i in range(df_peaks.shape[0]):
            for j in range(i+1, df_peaks.shape[0]):
                peak_table[i][j], neg_table[i][j] = self.modelScore(df_peaks, i, j)


        pt_size = len(peak_table) 
        # 2 ways to repeat values
        rownames = np.repeat(range(pt_size), pt_size)
        colnames = np.tile(range(pt_size), pt_size)
        # a mapping table of index and bin value (level)
        df_bintable = pd.DataFrame({'index': df_peaks.index, 'bin': df_peaks.bin})
        # a trick to flatten a 2d list
        values = sum(peak_table,[])
        neglects = sum(neg_table, [])
            
        # generate the dataframe
        df_models = pd.DataFrame({'base':rownames, 'thresh':colnames, 'score':values, 'neglects':neglects})

        # merge to get the bin value of base and thresh
        df_models = pd.merge(df_models, df_bintable, how='left', left_on='base', right_on='index')
        df_models = df_models.drop('index',1)
        df_models = df_models.rename(columns={'bin':'base_bin'})
        df_models = pd.merge(df_models, df_bintable, how='left', left_on='thresh', right_on='index')
        df_models = df_models.drop(['base', 'thresh','index'],1)
        df_models = df_models.rename(columns={'bin':'thresh_bin'})


        # get values of each bin and sum values of each model
        df_models = pd.merge(df_models, df_peaks[['value','bin']], how='left', left_on='base_bin', right_on='bin')
        df_models = df_models.drop('bin',1).rename(columns={'value':'base_value'})
        df_models = pd.merge(df_models, df_peaks[['value','bin']], how='left', left_on='thresh_bin', right_on='bin')
        df_models = df_models.drop('bin',1).rename(columns={'value':'thresh_value'})
        df_models['model_value'] = df_models.base_value + df_models.thresh_value

        # sort and ranking
        df_models = df_models.sort_values(by=['model_value','base_bin','thresh_bin'], ascending=False).sort_values(by='score')
        df_models['rank'] = range(pt_size**2)
        
        # print peak_table
        if showtable:
            print(pd.DataFrame(peak_table))
        
        # remove na from un-computed combinations in peak_table
        df_models = df_models.dropna()
        # mirror the values for the other half of the peak_table
        ##### Reason ####
        # Because one combination is computed only once, the 2 different sets of this combination ([a,b], [b,a]) all 
        # get values in different models. When later summing the values, they won't be able to be summed properly.
        # By mirroring the values, later summing can get a value for every combiniation in every model.
        # mirror = pd.DataFrame({'score':df_models.score, 'rank':df_models['rank'], 'base_bin':df_models.thresh_bin, 'thresh_bin':df_models.base_bin})
        mirror = df_models.rename(columns={'base_bin':'thresh_bin', 'thresh_bin':'base_bin', 'base_value':'thresh_value','thresh_value':'base_value'})

        mirrored_models = pd.concat([df_models, mirror], sort=False).sort_values(by='rank').reset_index(drop=True)
        levelsums= pd.DataFrame(mirrored_models.apply(lambda x: self.levelScore(df_peaks, x[1],x[2]), axis=1).tolist(), columns=['levelScore', 'dupLevels', 'delLevels'])
        mirrored_models = pd.concat([mirrored_models,levelsums], axis=1)

        # compute mean and std of levelscore
        # level_mean = np.mean(mirrored_models.levelScore)
        # level_std = np.std(mirrored_models.levelScore)
        # level_thresh = max(abs(level_mean - level_std), abs(level_mean + level_std))
        # if level_thresh < slcap:
        #     level_thresh = slcap
        # mirrored_models = mirrored_models[(abs(mirrored_models.levelScore) < level_thresh) ]

        # mirrored_models = mirrored_models[(mirrored_models.dupLevels<lcap) & (mirrored_models.delLevels>-lcap) & (abs(mirrored_models.levelScore) < slcap) ]        


        # combine and sort
        return mirrored_models






    # Find the closest number to n in ls
    def closestNum(self, n, ls):
        min_dif = 100
        idx = 0
        closest_value = 0
        for i in ls:
            dif = abs(n - i)
            if dif < min_dif:
                min_dif = dif
                idx = ls.index(i)
                closest_value = i
        return closest_value



    # Use 5 different interval settings to compute model scores, then combine 5 results to make final ranking.
    # The tricky step is to align bin(level) values, which are close but different in each setting.
    # Use the middle setting as the reference
    #
    #### Params ####
    # segments: the input, a list of dictionary. keys: chrom, start, end, probes, value
    # bins_per_interval: number of intervals of the middle setting, default=20
    # interval_step: the distance of two interval settings, default=4
    def integrateModels(self, segments):

        # interval settings except the bins_per_interval
        intervals = [self.bins_per_interval+self.interval_step, self.bins_per_interval+self.interval_step*2, 
            self.bins_per_interval+self.interval_step*3, self.bins_per_interval+self.interval_step*4]

        # model scores of bins_per_interval
        df_peaks_ref = self.computePeaks(segments)
        self.total_probes = sum(df_peaks_ref['value'])
        df_models = self.computeModels(df_peaks_ref)
        
        # combine results of other intervals to the central interval
        for i in intervals:
            df_interval_peaks = self.computePeaks(segments, bins_per_interval=i)

            if len(df_interval_peaks)<2:
                continue

            # a list of corresponding values of bins_per_interval in this interval
            # in the same order of bins_per_interval
            diffs = []
            # label names
            label_c = 'b'+str(self.bins_per_interval)
            # label_i = 'b'+str(i)
            # fill diffs
            for b in df_peaks_ref.bin:
                diffs.append(self.closestNum(b, list(df_interval_peaks.bin)))    
            # a mapping table of values in bins_per_interval and this interval
            bin_table = pd.DataFrame({label_c : df_peaks_ref.bin,                                   'base_bin' : diffs})
            
            # compute model scores of this interval and append the corresponding bin(level) value in the bins_per_interval
            # do for both base and thresh values

            # df_models_ext = self.computeModels(df_interval_peaks)

            df_models_ext = pd.merge(self.computeModels(df_interval_peaks), bin_table, how='left', on='base_bin')
            bin_table.rename(columns={'base_bin':'thresh_bin'}, inplace=True)
            df_models_ext = pd.merge(df_models_ext, bin_table, how='left', on='thresh_bin',                                       suffixes=('_bin', '_thresh'))    
            
            # merge info to the central interval
            df_models = pd.merge(df_models,                                         
                df_models_ext.loc[:,['score', 'rank',label_c+'_bin',label_c+'_thresh']].rename(columns={label_c+'_bin':'base_bin', label_c+'_thresh':'thresh_bin'}) ,                                                
                how='left', on=['base_bin', 'thresh_bin'],suffixes=('','_'+str(i)))
        
        # now can remove duplicate combiniations
        # df_models = df_models.drop_duplicates(subset='rank').reset_index(drop=True)

        # df_models['dist'] = abs(df_models.base_bin - df_models.thresh_bin)
        # if len(df_models[(df_models.dist >=self.thresh_min) & 
        #     (df_models.dist <= self.thresh_max)]) >0:
        #     df_models = df_models[(df_models.dist >=self.thresh_min) & 
        #     (df_models.dist <= self.thresh_max)]
        
        return df_models

        
    def integrateScores(self,segments):
        
        df_models = self.integrateModels(segments)


        # filter by levelScores
        level_mean = np.mean(df_models.levelScore)
        level_std = np.std(df_models.levelScore)
        level_thresh = max(abs(level_mean - level_std), abs(level_mean + level_std))
        if level_thresh < self.minlevel:
            level_thresh = self.minlevel
        df_models = df_models[(abs(df_models.levelScore) <= level_thresh) ]


        # filter by removing duplicate combinations
        # df_models = df_models.drop_duplicates(subset='rank').reset_index(drop=True)


        # filter by neglects ratio
        df_models['neglects_ratio'] = df_models.neglects / self.total_probes
        df_models = df_models[df_models.neglects_ratio < self.neglect_thresh]


        # filter by thresh distance
        df_models['dist'] = abs(df_models.base_bin - df_models.thresh_bin)
        if len(df_models[(df_models.dist >=self.thresh_min) & 
            (df_models.dist <= self.thresh_max)]) >0:
            df_models = df_models[(df_models.dist >=self.thresh_min) & 
            (df_models.dist <= self.thresh_max)]


            
        # filter by NA rows
        if df_models.isnull().values.any():
            if self.verbose:
                print('Remove NA rows:')
                print(df_models[df_models.isnull().any(axis=1)])
            df_models = df_models.dropna().reset_index(drop=True)

        
        # average rank
        df_models['ave_rank'] = df_models.filter(regex='rank*').sum(axis=1)/5

        # average score
        df = df_models.filter(regex='score*')
        df_models['ave_score'] = df.sum(axis=1)/5

        # average scaled score
        df -= df.min()
        df /= df.max()                                                                               
        df_models['ave_SS'] = df.sum(axis=1)/5
        
        # add ranking of each average
        df_models.sort_values(by='ave_rank', inplace=True)
        df_models['ar_rank'] = range(len(df))

        df_models.sort_values(by='ave_score', inplace=True)
        df_models['as_rank'] = range(len(df))

        df_models.sort_values(by='ave_SS', inplace=True)
        df_models['ass_rank'] = range(len(df))    

        if self.ranking_method == 'ar':
            df_models.sort_values(by='ave_rank', inplace=True)
        elif self.ranking_method == 'as':
            df_models.sort_values(by='ave_score', inplace=True)
         
        return df_models.reset_index(drop=True)






    # # determine the final baseline using consensus strategy
    # def determineBase(self, df_models):
    #     base_ranks = df_models.loc[df_models['ass_rank']<3,['base_bin', 'ass_rank']]
    #     thresh_ranks = df_models.loc[df_models['ass_rank']<3,['thresh_bin', 'ass_rank']]
    #     thresh_ranks = thresh_ranks.rename(columns={'thresh_bin':'base_bin'})
    #     bin_ranks = pd.concat([base_ranks,thresh_ranks])
    #     bin_ranks = bin_ranks.rename(columns={'base_bin':'bin'})
        
    #     # by score average, but if a bin only shows up once in rank 0, it always get highest score. Solution?
    # #     bin_ranks = bin_ranks.groupby('bin').mean().sort_values(by='ass_rank')
    #     # try counting, but when there are equal counts, what to to do? 
    # #    bin_ranks = bin_ranks.groupby('bin').count().sort_values(by='ass_rank', ascending=False)

    #     # combine both, comparison oder, count, value(probes)
    #     bin_ranks_mean = bin_ranks.groupby('bin').mean().rename(columns={'ass_rank':'mean'})
    #     bin_ranks_count = bin_ranks.groupby('bin').count().rename(columns={'ass_rank':'count'})
        
    #     return pd.merge(bin_ranks_mean, bin_ranks_count, how='left', on='bin')





    # determine the thresh level 
    def determineThresh(self, df_models):
        top_model = df_models.iloc[0,:]
        return round(abs(top_model['base_bin'] - top_model['thresh_bin']), 2)



    def determineBaseline(self, models):

        base_scores = {models.loc[0,'base_bin']:{}, models.loc[0,'thresh_bin']:{}}

        # levelScore
        if len(models) > 1 and models.loc[0,'rank'] == models.loc[1,'rank']:
            # compute level thresh
            level_mean = np.mean(models.levelScore)
            level_std = np.std(models.levelScore)
            level_thresh = max(abs(level_mean - level_std), abs(level_mean + level_std))
            if level_thresh < self.minlevel:
                level_thresh = self.minlevel

            ref_level = max(level_thresh - models.loc[0,'levelScore'], level_thresh - models.loc[1,'levelScore'] )
            base_scores[models.loc[0,'base_bin']]['level_score'] = (level_thresh - models.loc[0,'levelScore']) / ref_level
            base_scores[models.loc[1,'base_bin']]['level_score'] = (level_thresh - models.loc[1,'levelScore']) / ref_level
        else:
            base_scores[models.loc[0,'base_bin']]['level_score'] = 1
            base_scores[models.loc[0,'thresh_bin']]['level_score'] = 0
        


        # highest probe values
        ref_value = max(models.loc[0,'base_value'], models.loc[0,'thresh_value'])
        base_scores[models.loc[0,'base_bin']]['value_score'] = models.loc[0,'base_value'] / ref_value
        base_scores[models.loc[0,'thresh_bin']]['value_score'] = models.loc[0,'thresh_value'] / ref_value


        # closet to 2
        ref_dis = max(abs(2-abs(models.loc[0,'base_bin']-2)), abs(2-abs(models.loc[0,'thresh_bin']-2)))
        base_scores[models.loc[0,'base_bin']]['dis_score'] = abs(2-abs(models.loc[0,'base_bin']-2)) / ref_dis
        base_scores[models.loc[0,'thresh_bin']]['dis_score'] = abs(2-abs(models.loc[0,'thresh_bin']-2)) / ref_dis

        df_bscores = pd.DataFrame(base_scores)
        sums = df_bscores.sum(0)/3
        sums.name = 'sum'
        df_bscores = df_bscores.append(sums)

        baseline = sums[sums == max(sums)].index[0]

        return baseline, df_bscores

    # main
    def run(self, segments):
        peaks = self.computePeaks(segments, plot=self.plot)
        if len(peaks) >1:
            models = self.integrateScores(segments)
            models = models.round(6)
            if len(models) == 0:
                if self.verbose:
                    print('No models.')
                return ('No models.',)



            base_bin, best_vote = self.determineBaseline(models)


            thresh = self.determineThresh(models)
            
            levels = peaks.loc[:,['bin', 'value']]
            levels['R_level'] = round((levels['bin'] - base_bin) / thresh).astype(int)
            
            if self.verbose:
                print(peaks)
                print(models)
                # print('Base table:\n {}'.format(basedf))
                print("Base: {}".format(base_bin))
                print('Thresh: {}'.format(thresh))

            if self.outpath:
                with open(os.path.join(self.outpath, 'base_thresh.txt'), 'w') as fo:
                    # print('Interval with minimum Σe:\t{}'.format(base_candidates), file=fo)
                    print('Suggested baseline:\t{}'.format(base_bin), file=fo)
                    print('Suggested thresh:\t{}'.format(thresh), file=fo)

                models.to_csv(os.path.join(self.outpath, 'models.tsv'), sep='\t', index=False, float_format='%.4f')
                # basedf.to_csv(os.path.join(self.outpath, 'candidates.tsv'), sep='\t', index=False, float_format='%.2f')
                levels.to_csv(os.path.join(self.outpath, 'peaks.tsv'), sep='\t', index=False, float_format='%.2f')


            return (base_bin, thresh, models, best_vote, levels)
        else:
            if self.verbose:
                print('Not enough aberrant segments.')
            return ('Not enough aberrant segments.',)




