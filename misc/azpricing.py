#!/usr/bin/env python3

# Orginial developed by @bmoxon

""" Tool to query Azure pricing API for VM pricing information
args: region, family, version, 
"""
import argparse
import requests
import json, operator
import sys
import os
from tabulate import tabulate
import logging
import pandas as pd

logging.basicConfig(filename='azpricing.log', filemode='w', level=logging.DEBUG)

def parse_args():
    parser = argparse.ArgumentParser(description='Query Azure pricing API for VM pricing information')
    parser.add_argument('-r', '--region', help='Azure region to query, e.g. eastus', required=False)
    parser.add_argument('-f', '--family', help='VM family to query, e.g. ND, ND96asr', required=False)
    parser.add_argument('-v', '--version', help='VM family version, e.g. v4', required=False)
    parser.add_argument('-m', '--meter', help='Meter, e.g. "Spot" or "Low Priority"', required=False)
    return parser.parse_args()    

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

def addto_pricing_df(pricing_df, json_data):
    df = pd.DataFrame.from_dict(json_data['Items'])
    if df.empty:
        return pricing_df
    #drop Windows and Dedicated offerings
    trimmed_df = df[(df['productName'].str.contains('Windows') == False) & (df['productName'].str.contains('Dedicated') == False)]
    trimmed_df = trimmed_df.filter(items=['armSkuName', 'retailPrice', 'armRegionName', 'meterName', 'productName'])
    upd8_pricing_df = pd.concat([pricing_df, trimmed_df], ignore_index=True, sort=False)
    return upd8_pricing_df

def main():
    args = parse_args()
    region=args.region
    fam=args.family
    ver=args.version
    meter=args.meter
    
    pricing_df = pd.DataFrame(columns=['armSkuName', 'retailPrice', 'armRegionName', 'meterName', 'productName'])
    
    api_url = "https://prices.azure.com/api/retail/prices?api-version=2023-01-01-preview"
    query = construct_query(region, fam, ver, meter)
    response = requests.get(api_url, params={'$filter': query})

    if response == None:
        print("No response from API")
        sys.exit(1)
    elif response.status_code != 200:
        print("API returned status code: " + str(response.status_code))
        sys.exit(1)
        
    json_data = json.loads(response.text)
    final_pricing_df = addto_pricing_df(pricing_df, json_data)
    nextPage = json_data['NextPageLink']
    while(nextPage):
        response = requests.get(nextPage)
        json_data = json.loads(response.text)
        nextPage = json_data['NextPageLink']
        final_pricing_df = addto_pricing_df(final_pricing_df, json_data)

    final_pricing_df.sort_values(['armRegionName', 'retailPrice'], ascending=[1,1], inplace=True)
    print(tabulate(
        final_pricing_df.filter(items=['armSkuName', 'retailPrice', 'armRegionName', 'meterName']), \
        headers = 'keys', tablefmt = 'psql', showindex=False))
    
if __name__ == '__main__':
    main()