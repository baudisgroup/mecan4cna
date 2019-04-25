import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Return: mean, weighted mean and median
def medianValue(segment):
    seg_values = []
    for seg in segment:
        value = round(2**seg['value'] * 2, 4)
        probes = int(seg['probes'])
        seg_values.append([value, probes, np.log10(probes)])
    df = pd.DataFrame(seg_values, columns=['value','probes','log_probes'])

    percentiles = np.percentile(df.value, [20,80])

    df1 = df[(df.value > percentiles[0]) & (df.value < percentiles[1] )]

    return ( round(np.mean(df1.value) ,4),
             round(np.average(df1.value, weights=df1.log_probes), 4),
             round(np.median(df1.value), 4)) 

# plot all segments
def plotSegments(segments, filepath=None, probe_thresh=3):
    
    # do not plot segments less than probe_thresh
    
    # seperate segments into each chromosome
    chros = {}
    for seg in segments:
        chro = seg['chro']
        if chro not in chros:
            chros[chro] = [[int(seg['start']),int(seg['end']),int(seg['probes']),float(seg['value'])]]
        else:
            chros[chro].append([int(seg['start']),int(seg['end']),int(seg['probes']),float(seg['value'])])

    # determine the number of chromosomes
    num_plots = len(chros.keys())        
    plot_pos = 1
    plt.figure(figsize=(50,5))
    for k,v in chros.items():
        if plot_pos == 1:
            plt1 = plt.subplot(1, np.ceil(num_plots/1),plot_pos)
    #         plt1.set_yticks([-2,-1, 0, 1, 2])
        else:
            plt.subplot(1, np.ceil(num_plots/1),plot_pos, sharey=plt1)
    #         plt.yticks([])
        chro_len = 0
        for seg in v:
            if seg[2] > probe_thresh:
                x = [seg[0], seg[1]]
                yval = round(2**seg[3] * 2, 4)
                
                # supress extreme value
                if yval > 6:
                    yval = 6
                y = [yval, yval]
                
                # color of the line, hard coded
                if y[0] > 2.5:
                    color = 'b'
                elif y[0] < 1.5:
                    color = 'r'
                else:
                    color = 'k'
                plt.plot(x, y, color)
                plt.xticks([],[])
                chro_len = max(chro_len, seg[1])
    #     plt.plot([1,chro_len],[2,2], 'k:')
        plt.xlabel('chr'+str(k))
    #     plt.ylim([-2,2])
    #     plt.yticks([-2,-1, 0, 1, 2])
    #     plt1.set_yticks([-2,-1, 0, 1, 2])
        plot_pos +=1
    if filepath:
        plt.savefig(filepath, bbox_inches='tight')
    else:
        plt.show()
    

# plot a specific chromosome    
def plotChrom(segments, chrom, filepath=None, probe_thresh=3):
    
    # do not plot segments with less this number of probes
    
    chros = {}
    plt.figure(figsize=(30,5))
    for seg in segments:
        chro = seg['chro']
        if chro not in chros:
            chros[chro] = [[int(seg['start']),int(seg['end']),int(seg['probes']), round(2**seg['value'] * 2, 4),float(seg['value'])]]
        else:
            chros[chro].append([int(seg['start']),int(seg['end']),int(seg['probes']), round(2**seg['value'] * 2, 4),float(seg['value'])])
    for seg in chros[chrom]:
        if seg[2] > probe_thresh:
            x = [seg[0], seg[1]]
            yval = round(2**seg[4] * 2, 4)
            if yval > 6:
                yval = 6
            y = [yval, yval]
            if y[0] > 2.5:
                color = 'b'
            elif y[0] < 1.5:
                color = 'r'
            else:
                color = 'k'
            plt.plot(x, y, color)
    if filepath:
        plt.savefig(filepath, bbox_inches='tight')
    else:
        plt.show()
    return chros[chrom]

def file2list(fin):
    ls = []

    next(fin)
    for line in fin:
        line = line.split()
        ls.append({'chro':str(line[1]), 'start':int(line[2]), 'end':int(line[3]), 'probes':int(line[4]),'value':float(line[5])})
            
    return ls
            
