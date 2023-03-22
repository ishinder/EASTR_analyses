# get metadata from a filename: [sampleID][sep][treatment].[ext]
import os
import re

# remove path and extension from name
def parse_basename(filestr, extension):
    return re.split(extension, os.path.basename(filestr))[0]
# extracts the sampleID
# Precondition: samples are in the format "filepath/sampleID.{original,filtered}.bam"
def parse_sampleID(filestr=str, extension=str, delimiter=str):
    # extract filestr without the directory
    basename = parse_basename(filestr=filestr, extension=extension)
    sampleID = basename.split(delimiter)[0]
    # extract sampleID from filestr, assumed to be the format described by the header comments
    #sampleID = basename.split(file_sep)[0]
    return(sampleID)

# extracts the "treatment" ("original" or "filtered")
# Precondition: samples are in the format "filepath/sampleID.{original,filtered}.bam"
def parse_full_treatment(filestr=str, extension=str, delimiter=str):
    # extract filestr without the directory
    # extract sampleID from filestr, assumed to be the format described by the header comments
    basename = parse_basename(filestr=filestr, extension=extension)
    sampleID = parse_sampleID(filestr, extension=extension, delimiter=delimiter)
    treatment_start = len(sampleID) + len(delimiter)

    treatment = basename[treatment_start:] #os.path.basename(bam_filestr).split(file_sep)[1]
    return(treatment)


# assumes samples are in the format "filepath/sampleID.{original,filtered}.bam"
def infer_samples(file_list, t, delimiter, extension):
    samples = list()
    for filestr in file_list[t]:
        sampleID = parse_sampleID(filestr, extension=extension, delimiter=delimiter)
        if sampleID == None:
            raise Exception("naming convention of bam files should be" + \
                "[optional filepath/][sampleID].{original,filtered}.bam")
        samples.append(sampleID)
    return samples # assuming it is identical to samples_filtered