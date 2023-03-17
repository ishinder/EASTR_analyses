import argparse
import re
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Parse EASTR log file and output alignment stats")
    parser.add_argument("infile", type=str, help="Input log file path")
    parser.add_argument("outfile", type=str, help="Output TSV file path")
    return parser.parse_args()


def main(arglist=None):
    args = parse_args()
    aligners = {"hisat": [], "star": []}

    # Parse the log file
    stats = {}
    with open(args.infile, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Check if line contains a filename:
            match = re.match("brrrm! Vacuuming BAM file debris in: (.+)", line)
            if match:
                filename = match.group(1)
                if not filename:
                    raise ValueError(f"Empty filename in log file: {args.infile}")

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
                stats[filename]["Spliced Alignments"] = int(match.group(1))

            match = re.match("Alignments flagged for removal: (\d+)", line)
            if match:
                stats[filename]["Removed Alignments"] = int(match.group(1))

    # Create a DataFrame from the stats
    dfs = []
    for aligner, files in aligners.items():
        aligner_stats = {}
        for file in files:
            aligner_stats[file] = stats[file]
            aligner_stats[file]["SampleID"] = file.split("/")[-1].split("_", 1)[0]
        df = pd.DataFrame.from_dict(aligner_stats, orient="index")
        df.index.name = "Filename"
        df["Aligner"] = aligner
        dfs.append(df)
    df = pd.concat(dfs, axis=0, ignore_index=False)
    df = df.reset_index().set_index(["Aligner", "SampleID"])
    df = df.sort_index(level=["Aligner", "SampleID"])
    df = df.drop("Filename", axis=1)
    df["Percent Removed"] = df.apply(
        lambda row: row["Removed Alignments"] / row["Spliced Alignments"] * 100, axis=1
    )

    # Export to a TSV file
    df.to_csv(args.outfile, sep="\t")

def test_main():
    infile = "/ccb/salz2/shinder/projects/EASTR_tests2/lieber_sra/eastr_run.log"
    outfile = "/ccb/salz2/shinder/projects/EASTR_tests2/lieber_sra/eastr_run_stats.tsv"



if __name__ == "__main__":
    main()
