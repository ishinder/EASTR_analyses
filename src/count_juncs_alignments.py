

from collections import defaultdict
import csv
import glob
import os
import shlex
import subprocess
import sys
import tempfile
import argparse
import pandas as pd
import re


def parse_arguments():
    parser = argparse.ArgumentParser(description='Run EASTR analyses on reference junctions/alignments')
    parser.add_argument('--gtf_path', type=str, required=True, help='Path to GTF file')
    parser.add_argument('--ref_fa', type=str, required=True, help='Path to reference FASTA file')
    parser.add_argument('--bowtie2_index', type=str, required=True, help='Path to Bowtie2 index')
    parser.add_argument('--out_summary_file', type=str, help='Output summary file')
    parser.add_argument('--basedir', type=str, help='Base directory for input/output files')
    parser.add_argument('--metadata_file', type=str, help='Metadata file')
    parser.add_argument('--eastr_log_file', type=str, help='eastr log file output')
    return parser.parse_args()

def parse_eastr_log_file(eastr_log_file):
    aligners = {"hisat": [], "star": []}

    # Parse the log file
    stats = {}
    with open(eastr_log_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Check if line contains a filename:
            match = re.match("brrrm! Vacuuming BAM file debris in: (.+)", line)
            if match:
                filename = match.group(1)
                if not filename:
                    raise ValueError(f"Empty filename in log file: {eastr_log_file}")

                stats[filename] = {}
                aligner = None
                for a in aligners:
                    if a in filename:
                        aligner = a
                        break
                if aligner is None:
                    raise ValueError(f"Cannot determine aligner from filename: {filename}")
                aligners[aligner].append(filename)

            # Parse stats from the log file
            match = re.match("Spliced alignments identified: (\d+)", line)
            if match:
                stats[filename]["Total Spliced Alignments"] = int(match.group(1))

            match = re.match("Alignments flagged for removal: (\d+)", line)
            if match:
                stats[filename]["Removed Alignments"] = int(match.group(1))

    # Create a DataFrame from the stats
    dfs = []
    for aligner, files in aligners.items():
        aligner_stats = {}
        for file in files:
            aligner_stats[file] = stats[file]
            aligner_stats[file]["SRR ID"] = file.split("/")[-1].split("_", 1)[0]
        df = pd.DataFrame.from_dict(aligner_stats, orient="index")
        df.index.name = "Filename"
        df["Aligner"] = aligner
        dfs.append(df)
    df = pd.concat(dfs, axis=0, ignore_index=False)
    df = df.reset_index().set_index(["Aligner", "SRR ID"])
    df = df.sort_index(level=["Aligner", "SRR ID"])
    df = df.drop("Filename", axis=1)
    df["Percent Removed Alignments"] = df.apply(
        lambda row: row["Removed Alignments"] / row["Total Spliced Alignments"], axis=1
    )

    return df

def extract_splice_sites_gtf(gtf_path:str) -> dict:
    trans = {}
    gtf =open(gtf_path, "r")

    for line in gtf:
        line = line.strip()
        if line.startswith('#'):
            continue
        chrom, source, feature, start, end, score, \
                strand, frame, attributes = line.split('\t')

        start, end = int(start), int(end)

        if feature != 'exon':
            continue
        
        if start > end:
            raise Exception("Start of region can not be greater than end of region for:\n",line)

        values_dict = {}
        for attr in attributes.split(';'):
            if attr:
                attr, _, val = attr.strip().partition(' ')
                values_dict[attr] = val.strip('"')

        if 'transcript_id' not in values_dict:
            raise Exception("Exon does not contain transcript ID\n")

        transcript_id = values_dict['transcript_id']
        gene_id = values_dict['gene_id']

        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, gene_id, [[start, end]]]
        else:
            trans[transcript_id][3].append([start, end])


    for tran, [chrom, strand, gene_id, exons] in trans.items():
            exons.sort()


    junctions = defaultdict(dict)
    for tran, (chrom, strand, gene_id, exons) in trans.items():
        for i in range(1, len(exons)):
            if 'transcripts' not in junctions[(chrom, exons[i-1][1], exons[i][0]-1, strand)]:
                junctions[(chrom, exons[i-1][1], exons[i][0]-1, strand)]['transcripts'] = [gene_id, [tran]]
            else:
                junctions[(chrom, exons[i-1][1], exons[i][0]-1, strand)]['transcripts'][1].append(tran) #intron bed coordinates
    
    return junctions


def get_junctions_from_bed(bed_file):
    junctions = defaultdict(int)
    with open(bed_file, "r") as bed:
        for line in bed:
            if len(line.strip().split("\t")) < 6:
                print(line)
                
            chrom, start, end, name, score, strand = line.strip().split("\t")[0:6]
            start = int(start)
            end = int(end)
            try:
                score = int(score)
            except ValueError:
                score = 0
            junctions[(chrom, start, end, strand)] = score
    return junctions


def run_EASTR_on_gtf(gtf_path, bowtie2_index, ref_fa, ofile):
    #run EASTR on the gtf
    cmd = f"eastr --gtf {gtf_path} -r {ref_fa} -i {bowtie2_index} --out_removed_junctions {ofile}"
    process = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    if process.returncode != 0: #what is the return code if it fails? 
        print(f"Error running EASTR. Return code: {process.returncode}")
        print(f"stdout: {out.decode('utf-8')}")
        print(f"stderr: {err.decode('utf-8')}")
        sys.exit(1)
    
    removed_ref_junctions = get_junctions_from_bed(ofile)
    return removed_ref_junctions


def run_vacuum_on_removed_alns(bam_file, bed_file):
    dirname = os.path.dirname(bam_file)
    name = bam_file.split("/")[-1].split("_EASTR_filtered")[0]
    cmd = f"vacuum -V -r {dirname}/{name}_removed_ref_alignments.bam -o {dirname}/{name}_removed_nonref_alignments.bam {bam_file} {bed_file}"
    process = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    if process.returncode != 0:
        print(f"Error running vacuum. Return code: {process.returncode}")
        print(f"stdout: {out.decode('utf-8')}")
        print(f"stderr: {err.decode('utf-8')}")
        sys.exit(1)
    return out.decode()


def group_junctions(original_juncs_bed, removed_juncs_bed, ref_gtf_juncs, gtf_removed_juncs):
    original_juncs = get_junctions_from_bed(original_juncs_bed)
    removed_juncs = get_junctions_from_bed(removed_juncs_bed)

    ref_juncs = [j for j in original_juncs if j in ref_gtf_juncs]
    novel_juncs = [j for j in original_juncs if j not in ref_gtf_juncs]
    ref_removed_juncs = [j for j in removed_juncs if j in gtf_removed_juncs]
    nonref_removed_juncs = [j for j in removed_juncs if j not in gtf_removed_juncs]

    return ref_juncs, novel_juncs, ref_removed_juncs, nonref_removed_juncs


def load_metadata(metadata_file):
    sample_data = {}
    with open(metadata_file, "r") as metadata:
        reader = csv.reader(metadata)
        header = next(reader) # Get the header row
        
        for row in reader:
            sample_id = row[0]
            sample_data[sample_id] = {}
            for i in range(1, len(header)):
                sample_data[sample_id][header[i]] = row[i]
    
    return header,sample_data


def get_ref_junc_df(basedir,gtf_path, bowtie2_index, ref_fa, metadata_file):

    ofile = tempfile.NamedTemporaryFile(delete=False,suffix=".bed")
    ofile.close()
    ofile = ofile.name
    
    header,sample_data = load_metadata(metadata_file)
    ref_gtf_juncs = extract_splice_sites_gtf(gtf_path)
    gtf_removed_juncs = run_EASTR_on_gtf(gtf_path, bowtie2_index, ref_fa, ofile)

    for sample_name in sample_data:
        for software in ['hisat','star']:
            removed_alns = glob.glob(f"{basedir}/output/BAM/{sample_name}*{software}*_EASTR_filtered_removed_alignments.bam")[0]
            removed_juncs_bed = glob.glob(f"{basedir}/output/removed_junctions/{sample_name}*{software}*_removed_junctions.bed")[0]
            original_juncs_bed = glob.glob(f"{basedir}/output/original_junctions/{sample_name}*{software}*_original_junctions.bed")[0]
            ref_juncs, novel_juncs, ref_removed_juncs, nonref_removed_juncs = group_junctions(original_juncs_bed, removed_juncs_bed, ref_gtf_juncs, gtf_removed_juncs)
            sample_data[sample_name][f"{software}_ref_juncs"] = ref_juncs
            sample_data[sample_name][f"{software}_novel_juncs"] = novel_juncs
            sample_data[sample_name][f"{software}_ref_removed_juncs"] = ref_removed_juncs
            sample_data[sample_name][f"{software}_nonref_removed_juncs"] = nonref_removed_juncs
            out = run_vacuum_on_removed_alns(removed_alns, ofile)

            lines = out.split('\n')
            for line in lines:
                if 'Alignments flagged for removal:' in line:
                    num_removed = int(line.split(':')[-1].strip())
                    sample_data[sample_name][f"{software}_num_ref_alignments_removed"] = num_removed

    os.remove(ofile)

    return sample_data, header

def get_ref_junc_counts(basedir,gtf_path, bowtie2_index, ref_fa, metadata_file):

    sample_data, header = get_ref_junc_df(basedir,gtf_path, bowtie2_index, ref_fa, metadata_file)
    cols = ['Aligner'] + header + ['Reference Junctions', 
                            'Removed Reference Junctions', 'Novel junctions', 'Removed Novel Junctions', 
                            'Removed Reference Alignments']
    
    df = pd.DataFrame(columns=cols)
    for sample_name, data in sample_data.items():
        for software in ['hisat', 'star']:
            num_ref_juncs = len(sample_data[sample_name][f"{software}_ref_juncs"])
            num_novel_juncs = len(sample_data[sample_name][f"{software}_novel_juncs"])
            num_ref_removed_juncs = len(sample_data[sample_name][f"{software}_ref_removed_juncs"])
            num_nonref_removed_juncs = len(sample_data[sample_name][f"{software}_nonref_removed_juncs"])
            num_ref_removed_alns = sample_data[sample_name][f"{software}_num_ref_alignments_removed"] 
            
            # Create a dictionary for each row of data
            row_dict = {'Aligner': software, 'SRR ID': sample_name,
                        'Reference Junctions': num_ref_juncs, 'Removed Reference Junctions': num_ref_removed_juncs, 
                        'Novel junctions': num_novel_juncs, 'Removed Novel Junctions': num_nonref_removed_juncs, 
                        'Removed Reference Alignments': num_ref_removed_alns}
            
            for col in header:
                if col != 'SRR ID':
                    row_dict[col] = data[col]
            
            # Append the row to the dataframe
            df = df.append(row_dict, ignore_index=True)

    df = df.set_index(['Aligner', 'SRR ID'])
    df = df.sort_index()
    
    return df

def main():
    args = parse_arguments()
    gtf_path = args.gtf_path
    ref_fa = args.ref_fa
    bowtie2_index = args.bowtie2_index
    basedir = args.basedir
    metadata_file = args.metadata_file
    eastr_log_file = args.eastr_log_file
    out_summary_file = args.out_summary_file

    if basedir is None:
        basedir = os.getcwd()

    if metadata_file is None:
        metadata_file = f"{basedir}/metadata.csv"

    if eastr_log_file is None:
        eastr_log_file = f"{basedir}/eastr_run.log"

    if out_summary_file is None:
        out_summary_file = f"{basedir}/eastr_juncs_alns_summary.tsv"


    # gtf_path="/ccb/salz1/mpertea/stringtie/paper/hg38c_protein_and_lncRNA.gtf"
    # ref_fa = "/ccb/salz7-data/genomes/hg38/hg38p13.fa"
    # bowtie2_index = "/ccb/salz8-2/shinder/bt2_hg38_noPARs_index/hg38mod_noPARs"
    # basedir="/ccb/salz2/shinder/projects/EASTR_tests2/lieber_sra"
    # eastr_log_file = f"{basedir}/eastr_run.log"
    # metadata_file = f"{basedir}/metadata.csv"
    # out_summary_file = f"{basedir}/eastr_juncs_alns_summary.tsv"

    df1 = parse_eastr_log_file(eastr_log_file)
    df2 = get_ref_junc_counts(basedir, gtf_path, bowtie2_index, ref_fa, metadata_file)
    
    merged_df = pd.merge(df2, df1, left_index=True, right_index=True)
    merged_df['Removed Novel Alignments'] = merged_df['Removed Alignments'] - merged_df['Removed Reference Alignments']

    col_order = merged_df.columns[:-9].to_list() + ['Reference Junctions', 'Removed Reference Junctions', 'Novel junctions', 'Removed Novel Junctions','Total Spliced Alignments',
               'Removed Alignments', 'Removed Reference Alignments', 'Removed Novel Alignments', 'Percent Removed Alignments']
    
    merged_df[col_order].to_csv(out_summary_file,sep='\t')

if __name__ == '__main__':
    print("Starting script...")
    main()