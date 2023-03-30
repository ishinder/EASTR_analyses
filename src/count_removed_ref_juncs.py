from collections import defaultdict
import glob
import os
import argparse
import sys


def parse_arguments():
    parser = argparse.ArgumentParser(description='Run EASTR analyses on reference junctions/alignments')
    parser.add_argument('--gtf_path', type=str, required=True, help='Path to GTF file')
    parser.add_argument('--metadata_file', type=str, help='Metadata file')
    parser.add_argument('--basedir', type=str, help='Base directory for input/output files')
    parser.add_argument('--removed_juncs', type=str, help='Path to directory with removed junctions')
    parser.add_argument('--groupby', type=str, help='Group by column in metadata file')
    return parser.parse_args()


def find_removed_junctions_file(removed_juncs, aligner, sample_id):
    search_pattern = os.path.join(removed_juncs, f"{sample_id}*{aligner}*.bed")
    file_paths = glob.glob(search_pattern)
    if len(file_paths) == 0:
        raise ValueError(f"No removed junctions file found for sample {sample_id}")
    elif len(file_paths) > 1:
        raise ValueError(f"Multiple removed junctions files found for sample {sample_id}")
    return file_paths[0]


def main(metadata_file, basedir, removed_juncs, groupby, gtf_path):

    if basedir is None:
        basedir = os.getcwd()
    if removed_juncs is None:
        removed_juncs = os.path.join(basedir, 'output/removed_junctions')
    if metadata_file is None:
        metadata_file = f"{basedir}/metadata.csv"

    header, sample_data = load_metadata(metadata_file)
    if groupby:
        if groupby not in header:
            raise ValueError(f"Groupby column {groupby} not in metadata file")

    gtf_juncs = extract_splice_sites_gtf(gtf_path)
    all_bed_juncs = defaultdict(set) 
    groups = defaultdict(dict)
    for sample_id in sample_data:
        group = sample_data[sample_id][groupby] if groupby is not None else None
        
        for aligner in ['hisat', 'star']:
            bed_file = find_removed_junctions_file(removed_juncs, aligner, sample_id)
            juncs = get_junctions_from_bed(bed_file)
            all_bed_juncs[aligner].update(juncs)
            if group is not None:
                if group not in groups[aligner]:
                    groups[aligner][group] = set()
                groups[aligner][group].update(juncs)

    for aligner in ['hisat', 'star']:
        for group in groups[aligner]:
            juncs = groups[aligner][group]
            num_juncs = len(juncs.intersection(gtf_juncs))
            print(f"Number of {aligner} references junctions removed, grouped by {group}: {num_juncs}")

        num_juncs = len(all_bed_juncs[aligner].intersection(gtf_juncs))
        print(f"Number of {aligner} references junctions removed: {num_juncs}")

def test_main():
    sys.path.append('src') 
    from count_juncs_alignments import get_junctions_from_bed,extract_splice_sites_gtf, load_metadata
    
    gtf_path="/ccb/salz1/mpertea/stringtie/paper/hg38c_protein_and_lncRNA.gtf"
    basedir="/ccb/salz2/shinder/projects/EASTR_tests2/lieber_sra"
    metadata_file=f"{basedir}/metadata.csv"
    removed_juncs=f"{basedir}/output/removed_junctions"
    groupby="Library"
    main(metadata_file, basedir, removed_juncs, groupby, gtf_path)


if __name__ == '__main__':
    #import other scripts in same directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.append(current_dir) 
    from count_juncs_alignments import get_junctions_from_bed,extract_splice_sites_gtf, load_metadata

    print("Counting reference junctions...")
    args = parse_arguments()
    basedir = args.basedir
    gtf_path = args.gtf_path
    metadata_file = args.metadata_file
    removed_juncs = args.removed_juncs
    groupby = args.groupby
    main(metadata_file, basedir, removed_juncs, groupby, gtf_path)