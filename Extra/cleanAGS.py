#!/usr/bin/env python
# -*- coding:utf -*-
"""
#Created by Esteban.

Thu Nov 20, 2014

"""

import pandas as pd
import urllib.request as url
from time import sleep
import re

_ags = pd.read_csv("./AGS.csv", sep="\t", header=None, names=['ags', 'name'])
_clean_ags = _ags[_ags['ags'] >= 10000000000]
_clean_ags[_clean_ags['ags'] == 20000000000] = 2     # Hamburg
_clean_ags[_clean_ags['ags'] == 40110000000] = 4011  # Bremen
_clean_ags[_clean_ags['ags'] == 40120000000] = 4012  # Bremerhaven
_clean_ags[_clean_ags['ags'] == 110000000000] = 11   # Berlin
_clean_ags.to_csv("./AGS-Gemeinden.csv")
AGS = _clean_ags['ags'].tolist()

BASE_LINK = """https://ergebnisse.zensus2011.de/#dynTable:statUnit=PERSON;absRel=ANZAHL;ags={ags_id};agsAxis=X;yAxis={constrain};table"""
DOWNLOAD_LINK = """https://ergebnisse.zensus2011.de/auswertungsdb/download?csv=dynTable&tableHash=statUnit=PERSON;absRel=ANZAHL;ags={ags_id};agsAxis=X;yAxis={constrain}&locale=EN"""

CONSTRAINS = {
    "FAMSTND_KURZ": "Marital status",
    "ALTER_AF": "Age (eleven classes of years)",
    "HHGROESS_KLASS": "Size of private household"}

BASE_NAME = "../Data/Gemeinden/{}-{}.csv"


def reformat(i):
    i = str(i)
    if len(i) <= 11 and i != "11":  # Berlin does not have a leading 0
        i = "0" + i
    return(i)


def chunks(l, n):
    for i in range(0, len(l), n - 1):
        yield ",".join([reformat(i) for i in l[i:i + n]])


def getData(constrain):
    """gea the data for the given constrain."""

    dat_AGS = chunks(AGS, 100)
    for num, ags_c in enumerate(dat_AGS):
        to_download = DOWNLOAD_LINK.format(ags_id=ags_c, constrain=constrain)
        to_download = to_download.replace(" ", "")
        download_name = "../Data/Gemeinden/{}-{}.csv".format(
            constrain, num)

        url.urlretrieve(to_download, filename=download_name)

        sleep(1)  # be nice

    return(num)


def formatData(file_name, index):
    """Give the data the appropriate format."""

    data = pd.read_csv(
        file_name, sep=";", header=5, na_values="-",
        skip_footer=7, engine='python', index_col=index,
        encoding="latin-1")

    data = data.transpose()
    new_index = data.index
    new_index = new_index.map(lambda x: re.sub('[^0-9]', '', x))

    data.set_index(new_index, inplace=True)

    for col in data.columns:
        data[col] = data[col].map(lambda x: str(x).lstrip('(').rstrip(')'))

    return(data)


def cleanData(files_num, constrain):
    """clean and merge data files."""

    index = CONSTRAINS[constrain]

    first_file_name = BASE_NAME.format(constrain, 0)
    data = formatData(first_file_name, index)

    for f in range(1, files_num):
        file_name = BASE_NAME.format(constrain, f)
        this_data = formatData(file_name, index)
        data = data.append(this_data)

    data.to_csv(BASE_NAME.format(constrain, "all"))


def main():
    """ main function."""
    for constrain in CONSTRAINS:
        files_num = getData(constrain)
        cleanData(files_num, constrain)

if __name__ == "__main__":
    main()
