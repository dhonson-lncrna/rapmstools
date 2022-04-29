import numpy as np
import pandas as pd
from scipy.stats import f_oneway, tukey_hsd

def crapome_filter(df, 
                   spec='mouse'):
    '''Removes common background from RAP-MS data'''
    
    # List of some common background. Capitalize for human samples
    background = ('Krt','Dsp','Act','TUB','Tub','Lyz','Jup','nan')
    if spec == 'human':
        background = [i.upper() for i in background[:-1]]
        background.append('nan')
    else:
        pass
    
    # Filter out annotated contaminants
    df = df[~df['Protein IDs'].str.contains('CON')]
    
    # Filter out other known background
    for i in background:
        df = df[~df['Gene names'].str.contains(i)]
        
    return df

def import_rap(path,
               feature,
               targets,
               reps,
               spec):
    '''Imports RAP-MS data and filters contaminants'''
    
    # Make lists for the columns to keep
    keep = ['Protein IDs','Gene names']
    
    reps = np.array(np.arange(reps) + 1, dtype=str)
    samples = np.array(np.meshgrid(targets, reps)).T.reshape(-1,2)
    samples = ['_'.join(s) for s in samples]
    samples = [feature + ' ' + s for s in samples]
    
    keep = keep + samples
    
    # Import the data
    data = pd.read_csv(path, sep='\t', usecols=keep)
    data = data.astype({'Gene names':str})
    
    # Filter crapome
    data = crapome_filter(data, spec=spec)
    data = data.drop('Protein IDs',axis=1)
    
    # If column selection removed some samples, remove proteins that no longer show up
    data['Sum'] = np.sum(data.loc[:,keep[2:]], axis=1)
    data = data[data['Sum'] > 0]
    data = data.drop('Sum', axis=1)
    
    return data

def rap_anova(data,
              reps,
              thresh=0.05):
    '''Performs ANOVA for all factors in a dataset imported with import_rap'''
    # Set the index to Gene names
    data = data.set_index('Gene names')
    
    # Rename the columns for easier splitting
    data.columns = [i.split(' ')[-1] for i in data.columns]
    
    # Identify the individual RNAs and create a dictionary for indexing
    rnas = np.unique([i.split('_')[0] for i in data.columns])
    
    # It's easier to build the sample names from scratch than to pull them from the columns
    # to prevent key-value pairing errors
    reps = np.array(np.arange(reps) + 1, dtype=str)
    samples = []
    for r in rnas:
        samp_temp = np.array(np.meshgrid(r,reps)).T.reshape(-1,2)
        samples.append(['_'.join(i) for i in samp_temp])
    
        
    rnadict = dict(zip(rnas,samples))
    
    # Add a column to store p-values
    data['OneWay_Pval'] = np.zeros(len(data))
    
    # Perform ANOVA for each factor
    for i in data.index:
        subvals = [list(data.loc[i, rnadict[k]]) for k in rnas]
        data.loc[i,'OneWay_Pval'] = f_oneway(*subvals)[1]
        
    if thresh == False:
        return data
    
    else:
        return data[data['OneWay_Pval'] <= thresh]
    
def rap_tukey(anova,
              reps,
              fname,
              thresh=0.05):
    '''Performs Tukey test on all factors in a DataFrame output from rap_anova.
    Returns the log10-transformed p-value for RNAs that passed the threshold or 0 for those 
    that did not'''
    # Log transform threshold for later calculation
    thresh = -1*np.log10(thresh)
    
    # Identify the individual RNAs and create a dictionary for indexing
    anova = anova.drop('OneWay_Pval', axis=1)
    rnas = np.unique([i.split('_')[0] for i in anova.columns])
    
    # It's easier to build the sample names from scratch than to pull them from the columns
    # to prevent key-value pairing errors
    reps = np.array(np.arange(reps) + 1, dtype=str)
    samples = []
    for r in rnas:
        samp_temp = np.array(np.meshgrid(r,reps)).T.reshape(-1,2)
        samples.append(['_'.join(i) for i in samp_temp])
        
    rnadict = dict(zip(rnas,samples))
    
    # Make a new DataFrame to contain the values
    tukls = [np.zeros(len(anova)) for r in rnas] + [list(anova.index)]
    
    tukdf = pd.DataFrame(dict(zip(np.append(rnas, 'Gene names'), tukls)))
    tukdf = tukdf.set_index('Gene names')
    
    # Perform Tukey test for each factor
    for i in anova.index:
        subvals = [list(anova.loc[i,rnadict[k]]) for k in rnas]
        tukey = -1*np.log10(tukey_hsd(*subvals).pvalue)
        
        # Remove self correlation (p = 1)
        tukey = np.sort(tukey)[:,1:]
        
        # Calculate median p-values
        tmed = np.median(tukey, axis=1)
        
        # Determine which factors pass threshold in all pairwise tests
        tbool = tukey >= thresh
        tall = tbool.all(axis=1)
        
        # Update the dataframe
        row = np.zeros(len(tall))
        for j, v in enumerate(tall):
            if v == False:
                pass
            else:
                row[j] = tmed[j]
                
        tukdf.loc[i,:] = row
        
    tukdf.to_excel(fname)
        
    return tukdf

def reps_avg(df,
            reps):
    '''Calculates average values for replicates in RAP-MS DataFrame'''
    # New columns
    df = df.set_index('Gene names')
    cols = ['Gene names'] + list(np.unique([i.split('_')[0] for i in df.columns]))
    rows = [list(df.index)] + [np.zeros(len(df)) for i in cols[1:]]
    
    # Empty dataframe
    avg = pd.DataFrame(dict(zip(cols,rows)))
    
    # Groups
    groups = []
    for i in cols[1:]:
        groups.append([i + '_' + str(j) for j in np.arange(reps)+1])
    grdict = dict(zip(cols[1:],groups))
    
    # Averages
    for i in grdict.keys():
        avg.loc[:,i] = list(np.mean(df.filter(grdict[i]),axis=1))
        
    return avg
        
        
    
    
    