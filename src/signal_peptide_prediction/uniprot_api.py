"""Fetch UniProtKB search results one page at a time."""

import gzip
import json
import sys
import time

import requests


SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
SEARCH_FIELDS = (
    "accession,organism_name,lineage,protein_name,protein_existence,"
    "sequence,ft_signal,ft_transmem,cc_subcellular_location"
)
PAGE_SIZE = 500
RETRYABLE_STATUS_CODES = (429, 500, 502, 503, 504)


def read_page(session, url, params):
    for attempt in range(3):
        try:
            response = session.get(url, params=params, timeout=45)
            response.raise_for_status()
            payload = response.content
            # UniProt can send a gzip file without a Content-Encoding header.
            if payload.startswith(b"\x1f\x8b"):
                payload = gzip.decompress(payload)
            return json.loads(payload), response
        except (
            requests.HTTPError,
            requests.ConnectionError,
            requests.Timeout,
            requests.exceptions.ChunkedEncodingError,
            requests.exceptions.ContentDecodingError,
            EOFError,
        ) as error:
            if attempt == 2:
                raise

            delay = 2 ** attempt
            if isinstance(error, requests.HTTPError):
                if error.response.status_code not in RETRYABLE_STATUS_CODES:
                    raise
                retry_after = error.response.headers.get("Retry-After", "")
                if retry_after.isdecimal():
                    delay = max(delay, int(retry_after))
            print(f"Retrying UniProt page in {delay}s ({error})", file=sys.stderr)
            time.sleep(delay)


def uniprot_entries(query, max_pages, metadata):
    params = {
        "query": query,
        "format": "json",
        "fields": SEARCH_FIELDS,
        "size": PAGE_SIZE,
        "compressed": "true",
    }
    url = SEARCH_URL
    records_read = 0
    page_number = 0

    with requests.Session() as session:
        while url:
            data, response = read_page(session, url, params)
            release = response.headers["X-UniProt-Release"]
            total = int(response.headers["X-Total-Results"])
            if page_number == 0:
                metadata.update(
                    source=response.url,
                    query=query,
                    uniprot_release=release,
                    matching_records=total,
                )
            elif (release, total) != (metadata["uniprot_release"], metadata["matching_records"]):
                raise ValueError("UniProt's release or result count changed during collection.")

            page_number += 1
            entries = data["results"]
            records_read += len(entries)
            yield from entries

            if page_number == 1 or page_number % 10 == 0:
                print(f"Read {records_read:,} / {total:,} UniProt records", file=sys.stderr)

            # The next-page URL already includes the query and UniProt's cursor.
            url = response.links.get("next", {}).get("url")
            params = None
            if not url and records_read != total:
                raise ValueError(f"Expected {total:,} records, read {records_read:,}.")
            if page_number == max_pages:
                break

    metadata.update(pages=page_number, complete=records_read == total)
