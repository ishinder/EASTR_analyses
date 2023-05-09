#run gffcompare
import glob
import shlex
import subprocess

import sys

#from EASTR import utils
import glob
import argparse
import os
import re
import pandas as pd
import parse_file_metadata, juxtapose_outcomes


def parse_args(arglist):
    parser = argparse.ArgumentParser(
        prog="EASTR tests: individual GFF compare",
        description="Run gffcompare for original and filtered GTF samples")

    parser.add_argument(
        "--dir",
        help="base directory",
        default=os.getcwd())

    # Add an argument for the GTF directory
    parser.add_argument(
        "-gtf_dir",
        help="GTF directory",
        type=str,
        default=None)

    parser.add_argument(
        "-p",
        help="number of processes",
        type=int)

    parser.add_argument(
        "-r",
        help="reference GTF for sensitivity/precision calculation")

    parser.add_argument(
        "--outdir",
        help="base directory to output gffcompare results",
        default=None)

    parser.add_argument(
        "--filters",
        help="keywords to filter gtf stats file by, comma separated list",
        default=None)

    parser.add_argument(
        "--glob_by",
        help="regular expression string to select specific GTF files",
        default=None)
    
    parser.add_argument(
        "--delimiter",
        help="delimiter in file name separating sample name from treatment, etc.",
        default='_'
    )

    return parser.parse_args(arglist)

def make_dir(path):
    directory = os.path.join(path)
    os.makedirs(directory,exist_ok=True)

def get_filelist(directory, extension):
    filelist = []
    print("----------Checking all files in: ", os.path.join(directory, "*" + extension))
    for name in glob.glob(os.path.join(directory, "*" + extension)):
        filelist.append(name)
    return filelist 

def run_gffcompare(gtf_path, outdir, ref_gtf, name):
    make_dir(outdir)
    os.chdir(outdir)
    if not os.path.exists(f"{outdir}/{name}.stats"):
        cmd = f"gffcompare -r {ref_gtf} {gtf_path} -o {name}"
        print(cmd)
        subprocess.run(shlex.split(cmd), check=True)


def get_gffcompare_line(feature:str, filename:str):
    # grep "Transcript level" /ccb/salz2/shinder/projects/EASTR_tests/maize/gffcompare//*/*stats | awk -F'original/|filtered/' '{print $2}' | sort -nk1,1 | column -t
    cmd=f"grep \"{feature}\" {filename} | sed -E 's/[ \t]+/ /g' | sed -e 's/( /(/g'"

    return os.popen(cmd).read()

def get_gffcompare_stat(feature, filename, output, regex_str, flag_misses=True):
    featureStat = None
    search_featureStat = re.search(regex_str, output)
    if search_featureStat:
        featureStat = float(search_featureStat.group(1))
    elif flag_misses:
        print("Warning: Feature", feature, "not found for filename", filename)
        print("\tRegex:", regex_str)
    return featureStat


def get_sample_featureCount(feature:str, filename:str):
    #cmd=f"grep \"{filter_by}\" {outdir}/*/*.stats | awk -F'original/|filtered/' '{{print $2}}' | sort -nk1,1 | column -t"
    output = get_gffcompare_line(feature, filename)
    #print(output)

    # get raw count
    feature_stats = dict()
    featureCount = int(get_gffcompare_stat(feature=feature, filename=filename, \
        output=output, regex_str = feature + r'[ ]*:[ ]*(\d+)/?'))
    feature_stats["count"] = featureCount

    # get fraction, if it is a valid stat (not always)
    featurePercentage = get_gffcompare_stat(feature=feature, filename=filename, \
        output=output, regex_str = feature + r'[ ]*:[ ]*.*?[(]([0-9]+[.]?[0-9]+)[%][)]', flag_misses=False)
    feature_stats["percentage"] = featurePercentage

    # e.g.: GffCompare lists "novel exons" as count_a/count_b
    #       count_a is the actual count of "novel exons"
    #       count_b is the "total exons," which we extract as the "statRef" (reference point for the statistic)
    featureStatRef = get_gffcompare_stat(feature=feature, filename=filename, \
        output=output, regex_str = feature + r'[ ]*:[ ]*\d+/(\d+)', flag_misses=False)
    if featureStatRef:
        featureStatRef = int(featureStatRef)
    feature_stats["statRef"] = featureStatRef

    return feature_stats

def get_sample_loci_info(feature, filename):
    output = get_gffcompare_line(feature, filename)
    #print(output)

    # get raw count
    feature_stats = dict()
    featureLoci = get_gffcompare_stat(feature=feature, filename=filename, \
        output=output, regex_str = feature + r'[ ]*:[ ]*\d+[ ]*in[ ]*(\d+)[ ]*loci')
    feature_stats["loci"] = featureLoci

    return feature_stats

def get_sample_stats(filename:str, extension=".stats", delimiter='_'):
    # gffcompare file structure
    level_names=["Transcript level", "Exon level", "Intron level"]
    features_to_count=["Missed exons", "Novel introns", "Novel exons", "Novel loci", "Matching loci", \
        "Query mRNAs", "Reference mRNAs", "Matching transcripts"]
    loci_info=["Query mRNAs", "Reference mRNAs"]
    cols_report=["sampleID", "treatment", "Novel introns", "Novel exons", "Novel transcripts", \
               "Matching introns", "Matching exons", "Matching transcripts", \
                "Intron level sensitivity", "Intron level precision", "Exon level sensitivity", "Exon level precision", \
                "Transcript level sensitivity", "Transcript level precision"]

    stats_raw = dict()
    sampleID = parse_file_metadata.parse_sampleID(filename, extension=extension, delimiter=delimiter) # re.search('R\d+',filename).group()
    stats_raw["sampleID"] = sampleID
    treatment = parse_file_metadata.parse_full_treatment(filename, extension=extension, delimiter=delimiter)
    stats_raw["treatment"] = treatment
    
    # retrieve sensitivity and precision, for each level_name in level_names
    for line in open(filename):
        for level_name in level_names:
            #level_col = level_name.split(" level")[0]
            if re.search(level_name, line):
                level_info = re.findall('\d+.\d',line)
                stats_raw[level_name + " " + "sensitivity"] = float(level_info[0])
                stats_raw[level_name + " " + "precision"] = float(level_info[1])
    # retrieve counts for each feature (e.g.: "Missed exons", "Novel introns")
    for feature in features_to_count:
        feature_col = feature
        feature_stats = get_sample_featureCount(feature=feature, filename=filename)
        stats_raw[feature_col] = feature_stats["count"]
        stats_raw["percentage: " + feature_col] = feature_stats["percentage"]
        stats_raw["statRef: " + feature_col] = feature_stats["statRef"]
    for feature in loci_info:
        feature_col = "loci: "+ feature
        feature_stats = get_sample_loci_info(feature=feature, filename=filename)
        stats_raw[feature_col] = feature_stats["loci"]
        #print(feature, ":", end=" ")
        #print(feature_stats["count"], feature_stats["percentage"])
    
    # calculate and report essential stats (see cols_report)
    stats_report = dict.fromkeys(cols_report)
    #   novel count (except for transcripts) and sensitivity, precision
    for key in stats_report.keys():
        if key in stats_raw:
            stats_report[key] = stats_raw[key]
    
    #   reference-matching
    stats_report["Matching introns"] = stats_raw["statRef: Novel introns"] - stats_raw["Novel introns"]
    stats_report["Matching exons"] = stats_raw["statRef: Novel exons"] - stats_raw["Novel exons"]
    stats_report["Matching transcripts"] = stats_raw["Matching transcripts"]
    stats_report["Novel transcripts"] = stats_raw["Query mRNAs"] - stats_report["Matching transcripts"]

    # more precise calculation of transcript-level precision recall
    stats_check = dict()
    stats_check["Transcript level sensitivity"] = stats_report["Transcript level sensitivity"]
    stats_check["Transcript level precision"] = stats_report["Transcript level precision"]
    stats_report["Transcript level sensitivity"] = 100 * stats_report["Matching transcripts"]/(stats_raw["Reference mRNAs"])
    stats_report["Transcript level precision"] = 100 * stats_report["Matching transcripts"]/(stats_raw["Query mRNAs"])
    for key in stats_check.keys():
        if abs(stats_check[key] - stats_report[key]) > 0.1:
            print(f"-----WARNING: Check math for {key} (ID: {sampleID})-----")
            print(f"\tExpected: {stats_check[key]}")
            print(f"\tRecalculated: {stats_report[key]}")

    return stats_report #[sampleID, treatment_full, level_names, sensitivity, precision] 

def get_stats(directory=str, extension=str, delimiter=str):
    filelist = get_filelist(directory, extension)
    stats = pd.DataFrame()
    empty=True
    for filename in filelist:
        sample_stats_dict = get_sample_stats(filename=filename, extension=extension, delimiter=delimiter)
        sample_stats = pd.DataFrame(sample_stats_dict, index=[0])
        if empty:
            stats = sample_stats
        else:
            stats = pd.concat([stats, sample_stats], ignore_index=True)
        empty = False
    
    #stats = pd.DataFrame(stats, columns=["sampleID", "treatment", "category", "transcripts_sensitivity", "transcripts_precision"]).sort_values("sampleID").reset_index(drop=True)
    return stats

def print_gffcompare_results(filter_by, outdir, glob_by):
    cmd=f"grep \"{filter_by}\" {outdir}/*/*{glob_by}*.stats | awk -F'original/|filtered/' '{{print $2}}' | sort -nk1,1 | column -t"
    os.system(cmd)
    #result = subprocess.check_output(cmd, shell=True)
    #print(result)

def main(arglist=None):
    args = parse_args(arglist)
    p = args.p
    basedir = args.dir
    gtfpath = args.gtf_dir
    ref_gtf = args.r
    outdir = args.outdir
    glob_by = args.glob_by
    filters = args.filters
    delimiter = args.delimiter

    if gtfpath is None:
        gtfpath = f"{basedir}/GTF"
    
    if outdir is None:
        outdir = f"{basedir}/gffcompare"

    if glob_by is None:
        glob_by= "*"

    if filters is None:
        filters="Transcript level, Intron level, Missed exons, Novel introns, Novel exons, Novel loci, Matching transcripts"

    if p is None:
        p = 1

    #basedir="/ccb/salz2/shinder/projects/EASTR_tests2/lieber_sra"
    #outdir=f"{basedir}/gffcompare"
    #ref_gtf="/ccb/salz1/mpertea/stringtie/paper/hg38c_protein_and_lncRNA.gtf"

    gtf_list = {}
    for t in ["original","filtered"]:
        odir = f"{outdir}/{t}"

        print("GffCompare of", t, "to:", odir)
        make_dir(odir) # utils.make_dir(odir)

        os.chdir(odir)
        gtf_list[t] = glob.glob(f"{gtfpath}/{t}/*{glob_by}*.gtf")

    gtf_files = gtf_list["original"] + gtf_list["filtered"]
    run_gffcompare_parallel(gtf_files, outdir, ref_gtf, p=p)

    # filters=["Novel loci"]
    for filter_by in filters:
        sys.stdout.write(f"*-*-*-*-*-*-*-* Filtering by \"{filter_by}\" *-*-*-*-*-*-*-*\n")
        sys.stdout.flush()
        print_gffcompare_results(filter_by, outdir, glob_by)
        sys.stdout.write('\n')
        sys.stdout.flush()
    
    # report stats by sample and treatment
    stats_by_sample = get_stats(directory=os.path.join(outdir, "*"), extension=".stats", delimiter=delimiter)
    print(stats_by_sample)
    summary_outfile = os.path.join(outdir, "gffcompare_stats_summary.tsv")
    print("----------Writing stats summary to:", summary_outfile)
    stats_by_sample.to_csv(summary_outfile, index=False,sep='\t')
    # write stats as table
    table_outfile = os.path.join(outdir, "gffcompare_table.csv")
    juxtapose_outcomes.juxtapose_by_sample(stats_by_sample=stats_by_sample.copy(), outfilepath=table_outfile)
    print("----------Wrote table to:", table_outfile)

    return gtfpath

def test_main():
    arglist=["--dir", "/ccb/salz2/shinder/projects/EASTR_tests2/lieber_sra",
                   "-r", "/ccb/salz1/mpertea/stringtie/paper/hg38c_protein_and_lncRNA.gtf", 
                   "--filters", "Novel loci, Novel exons, Novel introns, Missed exons, Intron level, Transcript level, Matching transcripts", 
                   "--glob_by", "",
                   "-p", "12"]
    main(arglist)


if __name__ == '__main__':
    main()
