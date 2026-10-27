#!/usr/bin/env python
# coding: utf-8
from urllib.parse import urljoin, urlparse, urlsplit
import pandas as pd
import time
import os
import requests
from bs4 import BeautifulSoup


base = 'https://ftp.pride.ebi.ac.uk/pride/data/archive/'
URLs = [base]
depth = 4


def _is_within_base(url_: str, base_: str) -> bool:
    """Ensure we do not crawl outside the intended subtree."""
    u = urlparse(url_)
    b = urlparse(base_)
    return (u.scheme, u.netloc) == (b.scheme, b.netloc) and u.path.startswith(b.path)


for _ in range(depth):
    new_urls = []
    for url in URLs:
        if not _is_within_base(url, base):
            continue

        r = requests.get(url, timeout=60)        
        print('Error =', r.raise_for_status(), url)

        soup = BeautifulSoup(r.text, "html.parser")
        for a in soup.select("a[href]"):
            href = a.get("href")
            if not href.startswith(('/','?')):
                new_urls.append(urljoin(url,href))
    URLs = new_urls


data = pd.DataFrame(URLs, columns=['url'])
# data['pxd'] = data.url.apply(lambda x: os.path.split(x.strip('/'))[-1])
# data = data[data.pxd.str.startswith('PXD')].copy()
# data.set_index('pxd', inplace=True)
data.to_csv('PRIDE-URLs-files.csv', index=False, encoding='utf-8')