import argparse
import pandas as pd
import juxtapose_outcomes

def add_file(data_all:pd.DataFrame, filename:str, by:str, sep:str):
    data_new = pd.read_csv(filename, sep=sep)
    if len(data_all.columns) == 0:
        data_all = data_new
        return data_all
    if (set(data_new.columns) != set(data_all.columns)):
        raise ValueError("Files contain inconsistent columns/data")
    data_all = pd.concat([data_all, data_new])
    return data_all

def create_comparison(filelist:list, outfilepath:str, by:str, file_sep:str):
    data_all = pd.DataFrame()
    for filename in filelist:
        data_all = add_file(data_all, filename, by, file_sep)
    data_all = data_all.drop_duplicates()
    juxtapose_outcomes.juxtapose_by_sample(stats_by_sample=data_all.copy(), outfilepath=outfilepath)
    print()
    print(f"Written to {outfilepath}")
    return

def main():
    parser = argparse.ArgumentParser(description='Make a nested table comparing samples from different treatment groups')
    
    parser.add_argument('files', nargs='+', help='Space-delimited list of all data files')
    parser.add_argument('-o', default="comparison.tsv", help='summary file comparing two treatment groups', type=str)
    parser.add_argument('--by', default="sampleID", help='column name to group data by, default=\"sampleID\"', type=str)
    parser.add_argument('--sep', default="\t", help='file delimiter', type=str)

    args = parser.parse_args()

    create_comparison(args.files, args.o, by=args.by, file_sep=args.sep)

    return

if __name__ == '__main__':
    main()