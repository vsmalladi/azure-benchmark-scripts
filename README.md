# Sentieon DNAseq and DNAscope Benchmarking on Azure

This repository contains scripts used for benchmarking Sentieon DNAseq and DNAscope pipelines on Azure, on an Ubuntu Machine.

### Instance Setup

The script at `misc/instance_setup.sh` performs initial setup of the instance and download/installation of software packages used in the benchmark.

### Input datasets

Input datasets for the benchmark are recorded in the `config/config.yaml` (and `config/config_arm.yaml`) files. With the exception of the Element dataset, that you will have to download on your own.

We recomend downloading all the files and placing them in azure blob storage. You can use [AzCopy](https://docs.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-blobs#copy-a-container-to-another-storage-account) to transfer the required files to your own Storage account [using a shared access signature](https://docs.microsoft.com/en-us/azure/storage/common/storage-sas-overview) with "Write" access. Then we recomend updated tine configs to use a shared access signature to each file. The pipeline will automatic download input files. Visiting the Element web page recorded in the config file using an internet browser is required for downloading the Element dataset.

### Running benchmarks

The script at `misc/run_benchmarks.sh` was used to run the benchmarks. Before running the scripts, the `SENTIEON_LICENSE` variable inside the script should be properly set.


### Collecting runtimes

The script at `misc/collect_runtimes.py` was used to collect bechmarking results. To get pricing use `misc/azpricing.py` which published information on the Azure pricing APIs.


[Azure Retail Prices overview](https://learn.microsoft.com/en-us/rest/api/cost-management/retail-prices/azure-retail-prices)

```
$ ./azpricing.py --help
usage: azpricing.py [-h] [-r REGION] [-f FAMILY] [-v VERSION] [-m METER]

Query Azure pricing API for VM pricing information

optional arguments:
  -h, --help            show this help message and exit
  -r REGION, --region REGION
                        Azure region to query, e.g. eastus
  -f FAMILY, --family FAMILY
                        VM family to query, e.g. ND, ND96asr
  -v VERSION, --version VERSION
                        VM family version, e.g. v4
  -m METER, --meter METER
                        Meter, e.g. "Spot" or "Low Priority"

$ ./azpricing.py -r eastus2 -f ND -v 4
+---------------------------+---------------+-----------------+-------------------------------+
| armSkuName                |   retailPrice | armRegionName   | meterName                     |
|---------------------------+---------------+-----------------+-------------------------------|
| Standard_ND96amsr_A100_v4 |        6.554  | eastus2         | ND96amsr A100 v4 Low Priority |
| Standard_ND96amsr_A100_v4 |       14.4188 | eastus2         | ND96amsr A100 v4 Spot         |
| Standard_ND96amsr_A100_v4 |       32.77   | eastus2         | ND96amsr A100 v4              |
+---------------------------+---------------+-----------------+-------------------------------+
```


