#run gffcompare
import glob
import shlex
import subprocess
from collections import defaultdict
import sys
from EASTR import utils
import glob
import argparse
import os

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

    return parser.parse_args(arglist)

def make_dir(path):
    directory = os.path.join(path)
    os.makedirs(directory,exist_ok=True)

def run_gffcompare(gtf_path, outdir, ref_gtf, name):
    make_dir(outdir)
    os.chdir(outdir)
    if not os.path.exists(f"{outdir}/{name}.stats"):
        cmd = f"gffcompare -r {ref_gtf} {gtf_path} -o {name}"
        subprocess.run(shlex.split(cmd), check=True)


def print_gffcompare_results(filter_by, outdir, glob_by):
    cmd=f"grep \"{filter_by}\" {outdir}/*/*{glob_by}*stats | awk -F'original/|filtered/' '{{print $2}}' | sort -nk1,1 | column -t"
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

    if gtfpath is None:
        gtfpath = f"{basedir}/GTF"
    
    if outdir is None:
        outdir = f"{basedir}/gffcompare"

    if glob_by is None:
        glob_by= "*"

    if filters is None:
        filters="Transcript level, Intron level, Missed exons, Novel introns, Novel exons, Novel loci, Matching transcripts"

    filters = filters.split(', ')

    basedir="/ccb/salz2/shinder/projects/EASTR_tests2/lieber_sra"
    outdir=f"{basedir}/gffcompare"
    ref_gtf="/ccb/salz1/mpertea/stringtie/paper/hg38c_protein_and_lncRNA.gtf"

    gtf_list = {}
    for t in ["original","filtered"]:
        odir = f"{outdir}/{t}"
        utils.make_dir(odir)
        os.chdir(odir)
        gtf_list[t] = glob.glob(f"{gtfpath}/{t}/*{glob_by}*.gtf")
        for file in gtf_list[t]:
            name=os.path.basename(file).split('.')[0]
            if not os.path.exists(f"{odir}/{name}.stats"):
                run_gffcompare(file, odir, ref_gtf, name)

    # filters=["Novel loci"]
    for filter_by in filters:
        sys.stdout.write(f"*-*-*-*-*-*-*-* Filtering by \"{filter_by}\" *-*-*-*-*-*-*-*\n")
        sys.stdout.flush()
        print_gffcompare_results(filter_by, outdir, glob_by)
        sys.stdout.write('\n')
        sys.stdout.flush()
    
    return gtfpath

def test_main():
    arglist=["--dir", "/ccb/salz2/shinder/projects/EASTR_tests2/lieber_sra",
                   "-r", "/ccb/salz1/mpertea/stringtie/paper/hg38c_protein_and_lncRNA.gtf", 
                   "--filters", "Novel loci, Novel exons, Novel introns, Missed exons, Intron level, Transcript level, Matching transcripts", 
                   "--glob_by", ""]


if __name__ == '__main__':
    main()
