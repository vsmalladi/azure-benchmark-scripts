#!/usr/bin/env python3

"""
Collect the DNAseq runtimes from one or more benchmark runs
"""

import argparse
import collections
import glob
import os.path
import sys

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input_directory", help="A directory containing benchmark results")
    parser.add_argument("--instance_price", type=float, help="The instance price in $/hr")
    parser.add_argument("--instance_type", help="The instance type")
    parser.add_argument("--outfile", default=sys.stdout, type=argparse.FileType('w'))
    return parser.parse_args(argv)

def bm_file_to_sec(bm_path):
    """ Extract the runtime in seconds from the snakemake benchmark file """
    fh = open(bm_path)
    _hdr = fh.readline()
    values = fh.readline().rstrip().split('\t')
    return float(values[0])

def main(args):

    targets=[
        "DNAseq/Illumina_HG002_HiSeqX_40x",
        "DNAseq/Illumina_HG002_NovaSeq_40x",
        "DNAscope/Illumina_HG002_HiSeqX_40x",
        "DNAscope/Illumina_HG002_NovaSeq_40x",
        "DNAscope/Element_HG002_100x",
        "DNAscope/Ultima_HG002_cram",
        "DNAscope_LongRead/PacBio_HG002_HiFi_Chem2",
    ]


    # Iterate over each target
    overall_results = collections.OrderedDict()
    in_dir = str(args.input_directory)
    price = args.instance_price
    instance_type = args.instance_type
    for target in targets:
        runtime_in_sec = collections.OrderedDict()
        price = float(price)

        target_split = target.split('/')
        pipeline = target_split[0] 
        sample = target_split[1]
        if (pipeline == "DNAseq"):

            # Alignment runtimes
            alignment_results = glob.glob(os.path.join(in_dir, "ramdisk/bwa_mem", sample,  "0/*/*/alignments.bam.benchmark.txt"))
            runtime_in_sec["Alignment"] = max([bm_file_to_sec(x) for x in alignment_results])  # One or more alignemnt jobs run in parallel. All jobs need to finish before the next step, so total runtime is the max of all runtimes

            # Preprocessing
            runtime_in_sec["LocusCollector"] = bm_file_to_sec(os.path.join(in_dir, "ramdisk/dedup/lc", sample, "score.txt.gz.benchmark.txt"))
            runtime_in_sec["Dedup"] = bm_file_to_sec(os.path.join(in_dir, "ramdisk/dedup/dedup", sample, "dedup.bam.benchmark.txt"))
            runtime_in_sec["QualCal"] = bm_file_to_sec(os.path.join(in_dir, "qualcal", sample, "recal.table.benchmark.txt"))
            
            # Variant Calling
            runtime_in_sec["Haplotyper"] = bm_file_to_sec(os.path.join(in_dir, "haplotyper", sample, "calls.vcf.gz.benchmark.txt"))
        
        elif (pipeline == 'DNAscope'):
            if (sample != 'Ultima_HG002_cram'):
                # Alignment runtimes
                alignment_results = glob.glob(os.path.join(in_dir, "ramdisk/bwa_mem", sample,  "0/*/*/alignments.bam.benchmark.txt"))
                runtime_in_sec["Alignment"] = max([bm_file_to_sec(x) for x in alignment_results])  # One or more alignemnt jobs run in parallel. All jobs need to finish before the next step, so total runtime is the max of all runtimes
            

                # Preprocessing
                runtime_in_sec["LocusCollector"] = bm_file_to_sec(os.path.join(in_dir, "ramdisk/dedup/lc", sample, "score.txt.gz.benchmark.txt"))
                runtime_in_sec["Dedup"] = bm_file_to_sec(os.path.join(in_dir, "ramdisk/dedup/dedup", sample, "dedup.bam.benchmark.txt"))
            
            # Variant Calling
            runtime_in_sec["Dnascope"] = bm_file_to_sec(os.path.join(in_dir, "dnascope", sample, "calls.vcf.gz.benchmark.txt"))
            runtime_in_sec["Dnascope_tmp"] = bm_file_to_sec(os.path.join(in_dir, "dnascope_tmp", sample, "calls.vcf.gz.benchmark.txt"))
        
        elif (pipeline == 'DNAscope_LongRead'):
            # Alignment runtimes
            alignment_results = glob.glob(os.path.join(in_dir, "ramdisk/minimap2", sample, "*/alignments.bam.benchmark.txt"))
            runtime_in_sec["Alignment"] = sum([bm_file_to_sec(x) for x in alignment_results])  # One or more alignemnt jobs run in parallel. All jobs need to finish before the next step, so total runtime is the max of all runtimes

            # Preprocssing
            runtime_in_sec["MergePb"] = bm_file_to_sec(os.path.join(in_dir, "ramdisk/merge_pb", sample, "merged.bam.benchmark.txt"))
            
            # Variant Calling
            runtime_in_sec["Dnascope"] = bm_file_to_sec(os.path.join(in_dir, "dnascope_lr", sample, "calls.vcf.gz.benchmark.txt"))
        
        # summary results
        results = collections.OrderedDict()
        results["Total (s)"] = sum(runtime_in_sec.values())
        results["Total (min)"] = sum(runtime_in_sec.values()) / 60
        results["$/hr"] = price
        results["Total compute cost"] = sum(runtime_in_sec.values()) / 60 / 60 * price

        # Collect everything together
        runtime_in_sec.update(results)
        overall_results[target] = runtime_in_sec.copy()
        


    # Output a nice table
    print(args.instance_type, file=args.outfile)
    print("\n", file=args.outfile)
    for sample in targets:
        row_names = list(overall_results[sample].keys())
        print(sample, file=args.outfile)
        for row_name in row_names:
            row = [row_name] + [str(overall_results[sample][row_name])]
            print('\t'.join(row), file=args.outfile)
        print("\n", file=args.outfile)

    return 0

if __name__ == "__main__":
    args = parse_args()
    main(args)
