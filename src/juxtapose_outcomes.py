import pandas as pd

def display_elements(iterable, sep=" , "):
    iterStr = ""
    for element in iterable:
        iterStr += element + sep
    iterStr = iterStr[:-len(sep)]
    return iterStr

# given conditions A and B and a single stats dataframe containing pairwise A,B data,
#   reorganize the data into a nested table juxtaposing A,B by individual (matching ID) for each statistic (dataframe column)
def juxtapose_by_sample(stats_by_sample:pd.DataFrame, outfilepath:str, pairedID_col="sampleID", treatment_col="treatment"):
    subcols = 3
    subcol_sep = ',,,'
    print("---------Juxtaposing stats")
    treatments = stats_by_sample[treatment_col].unique()
    treatments_tostr = display_elements(treatments)

    # user will explicitly specify what type of relative stats are required
    treatment_a = None
    treatment_b = None
    redo = 0
    while (treatment_a not in treatments) or (treatment_b not in treatments):
        if redo > 0:
            print("Entered invalid treatment(s). Please use listed options!")
        treatment_a = str(input("Choose reference/\"original\" condition from the following: " + treatments_tostr + "\n"))
        treatment_b = str(input("Choose manipulated or \"filtered\" condition from the following: " + treatments_tostr + "\n"))
        redo += 1
    print()

    with open(outfilepath, 'w') as outfile:
        header_cols=list()
        sub_cols = list()
        for col in stats_by_sample.columns:
            header_cols.append(col)
            sub_cols.append("control")
            sub_cols.append("modified")
            sub_cols.append("relative change (%)")
        outfile.write(subcol_sep.join(header_cols) + '\n')
        outfile.write(','.join(sub_cols) + '\n')
        IDs = stats_by_sample.loc[stats_by_sample[treatment_col] == treatment_a, pairedID_col].values
        for ID in IDs:
            # extract numerical data (ignore ID and treatment_col, which are known)
            ID_stats_a = stats_by_sample.loc[(stats_by_sample[pairedID_col] == ID) & (stats_by_sample[treatment_col] == treatment_a), \
                ~stats_by_sample.columns.isin([pairedID_col, treatment_col])].copy()
            ID_stats_b = stats_by_sample.loc[(stats_by_sample[pairedID_col] == ID) & (stats_by_sample[treatment_col] == treatment_b), \
                ~stats_by_sample.columns.isin([pairedID_col, treatment_col])].copy()

            # check valid experimental design (for this stats comparison, there should not be duplicate IDs per treatment)
            if len(ID_stats_a) > 1:
                raise Exception("Multiple entries found for " + ID + ", under treatment " + treatment_a +\
                ". Unless these results are identical, please differentiate or remove one")
            if len(ID_stats_b) > 1:
                raise Exception("Multiple entries found for " + ID + ", under treatment " + treatment_b +\
                ". Unless these results are identical, please differentiate or remove one")

            # juxtapose
            sample_data = list()
            sample_data.append(ID)
            for i in range(1, subcols):
                sample_data.append("")
            sample_data.append(treatment_a)
            sample_data.append(treatment_b)
            sample_data.append("relative_change")
            for col in ID_stats_a.columns: # must be ordered
                value_a = ID_stats_a[col].values[0]
                value_b = ID_stats_b[col].values[0]
                sample_data.append(str(value_a))
                sample_data.append(str(value_b))
                sample_data.append(str( round(100*(value_b - value_a)/value_a, 2) ))
            outfile.write(','.join(sample_data) + '\n')
    return