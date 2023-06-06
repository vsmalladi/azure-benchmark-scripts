#!/usr/bin/env python3

"""
Collect the DNAseq runtimes from one or more benchmark runs
"""

import argparse
import collections
import glob
import os.path
import sys
import requests
import json, operator
import sys
import os
import logging
import pandas as pd

pipelines = {
    "DNAseq": collections.OrderedDict([
        (
            "Alignment",
            {
                "path": "ramdisk/bwa_mem/{sample}/0/*/*/alignments.bam.benchmark.txt",
                "reduce_func": max,
            },
        ),
        (
            "LocusCollector",
            {
                "path": "ramdisk/dedup/lc/{sample}/score.txt.gz.benchmark.txt",
            },
        ),
        (
            "Dedup",
            {
                "path": "ramdisk/dedup/dedup/{sample}/dedup.bam.benchmark.txt",
            },
        ),
        (
            "QualCal",
            {
                "path": "qualcal/{sample}/recal.table.benchmark.txt",
            },
        ),
        (
            "Haplotyper",
            {
                "path": "haplotyper/{sample}/calls.vcf.gz.benchmark.txt",
            }
        ),
    ]),
    "DNAscope CRAM": collections.OrderedDict([
        (
            "DNAscope_tmp",
            {
                "path": "dnascope_tmp/{sample}/calls.vcf.gz.benchmark.txt",
            },
        ),
        (
            "DNAscope",
            {
                "path": "dnascope/{sample}/calls.vcf.gz.benchmark.txt",
            },
        ),
    ]),
    "DNAscope": collections.OrderedDict([
        (
            "Alignment",
            {
                "path": "ramdisk/bwa_mem/{sample}/0/*/*/alignments.bam.benchmark.txt",
                "reduce_func": max,
            },
        ),
        (
            "LocusCollector",
            {
                "path": "ramdisk/dedup/lc/{sample}/score.txt.gz.benchmark.txt",
            },
        ),
        (
            "Dedup",
            {
                "path": "ramdisk/dedup/dedup/{sample}/dedup.bam.benchmark.txt",
            },
        ),
        (
            "DNAscope_tmp",
            {
                "path": "dnascope_tmp/{sample}/calls.vcf.gz.benchmark.txt",
            },
        ),
        (
            "DNAscope",
            {
                "path": "dnascope/{sample}/calls.vcf.gz.benchmark.txt",
            },
        ),
    ]),
    "DNAscope LongRead": collections.OrderedDict([
        (
            "Alignment",
            {
                "path": "ramdisk/minimap2/{sample}/*/alignments.bam.benchmark.txt",
            },
        ),
        (
            "Merge",
            {
                "path": "ramdisk/merge_pb/{sample}/merged.bam.benchmark.txt",
            },
        ),
        (
            "DNAscope-LR",
            {
                "path": "dnascope_lr/{sample}/calls.vcf.gz.benchmark.txt",
            },
        ),
    ])
}

pipeline_samples = {
    "DNAseq": [
        "Illumina_HG002_HiSeqX_40x",
        "Illumina_HG002_NovaSeq_40x",
    ],
    "DNAscope": [
        "Illumina_HG002_HiSeqX_40x",
        "Illumina_HG002_NovaSeq_40x",
        "Element_HG002_100x",
    ],
    "DNAscope CRAM": [
        "Ultima_HG002_cram",
    ],
    "DNAscope LongRead": [
        "PacBio_HG002_HiFi_Chem2",
    ],
}

row_order = [
    "Sample", "Pipeline","Alignment", "Merge", "LocusCollector", "Dedup", "QualCal", "Haplotyper",
    "DNAscope_tmp", "DNAscope", "DNAscope-LR", "Total (s)", "Total (min)",
    "$/hr", "Total compute cost"
]

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input_directory", help="A directory containing benchmark results")
    parser.add_argument('-f', '--family', help='VM family to query, e.g. ND, ND96asr')
    parser.add_argument('-v', '--version', help='VM family version, e.g. v4')
    parser.add_argument('-r', '--region', help='Azure region to query, e.g. eastus')
    parser.add_argument('-m', '--meter', help='Meter, e.g. "Spot" or "Low Priority"')
    parser.add_argument("--outfile", help='CSV format file' , type=argparse.FileType('w'))
    return parser.parse_args(argv)

def construct_query(region, fam, ver, meter):
    # only VMs under consumption plans (PAYG, Low Priority, Spot); no reserved instances
    query = 'serviceName eq \'Virtual Machines\' and priceType ne \'Reservation\''
    if region != None:
        query += ' and armRegionName eq \'' + region + '\''
    if fam != None:
        query = query + ' and contains(skuName, \'' + fam + '\')'
    if ver != None:
        query = query + ' and contains(skuName, \'' + ver + '\')'
    if meter != None:
        query = query + ' and contains(meterName, \'' + meter + '\')'

    logging.debug('query = {}'.format(query))
    return query

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


    # Azure pricing
    region=args.region
    fam=args.family
    ver=args.version
    meter=args.meter
    
    api_url = "https://prices.azure.com/api/retail/prices?api-version=2023-01-01-preview"
    query = construct_query(region, fam, ver, meter)
    response = requests.get(api_url, params={'$filter': query})
    json_data = json.loads(response.text)
    df = pd.DataFrame.from_dict(json_data['Items'])
    #drop Windows and Dedicated offerings
    trimmed_df = df[(df['productName'].str.contains('Windows') == False) & (df['productName'].str.contains('Dedicated') == False)]
    instance_price = trimmed_df['retailPrice'].iloc[0]
    instance_type = trimmed_df['armSkuName'].iloc[0]

    if response == None:
        print("No response from API")
        sys.exit(1)
    elif response.status_code != 200:
        print("API returned status code: " + str(response.status_code))
        sys.exit(1)



    # Iterate over each target
    overall_results = pd.DataFrame()
    in_dir = str(args.input_directory)
    price = instance_price
    instance_type = instance_type
    for target in targets:
        runtime_in_sec = collections.OrderedDict()
        price = float(price)
        for pipeline, samples in pipeline_samples.items():
            for sample in samples:
                runtime_in_sec = collections.OrderedDict()
                missing_stage = False
                for stage, d in pipelines[pipeline].items():
                    path = d["path"].format(sample=sample)
                    fullpath = os.path.join(in_dir, path)
                    stage_results = glob.glob(fullpath)
                    if not stage_results:
                        missing_stage = True
                        print(f"Missing files for instance '{instance_type}', pipeline '{pipeline}', and sample '{sample}' at: {fullpath}", file=sys.stderr, flush=True)
                        break

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
        results['Sample'] = sample
        results['Pipeline'] = pipeline

        # Collect everything together
        runtime_in_sec.update(results)
        df_results = pd.DataFrame(runtime_in_sec, index=[0])
        overall_results = pd.concat([overall_results, df_results])
        


    # Output a nice table
    overall_results = overall_results[row_order]
    overall_results.to_csv(args.outfile, sep=',', index=False)


    return 0

if __name__ == "__main__":
    args = parse_args()
    main(args)
