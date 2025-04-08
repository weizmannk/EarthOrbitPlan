# coding: utf-8
"""
---------------------------------------------------------------------------------------------------
                                    ABOUT
@author         : Ramodgwendé Weizmann KIENDREBEOGO
@email          : kiend.weizman7@gmail.com / weizmann.kiendrebeogo@oca.eu
@repo           : https://github.com/weizmannk/ultrasat-gw-followup.git
@createdOn      : December 2023
@description    : A tool for interacting with the Zenodo API, facilitating the download of files based on DOI.
                This class provides functionality to retrieve the latest version DOI associated with a provided
                permanent DOI, and subsequently download the corresponding file from Zenodo.

                You can easily download another dataset from Zenodo by replacing the permanent_doi
                  with a new one.

Usage:
    python3 zenodo_downloader.py
---------------------------------------------------------------------------------------------------
"""

import os

import requests
from tqdm.auto import tqdm

# If the link to the data is required manually, it can be accessed at : https://doi.org/10.5281/zenodo.14142969

permanent_doi = "14142969"
file_name = "runs_SNR-10.zip"


class ZenodoDownloader:
    """A class to interact with the Zenodo API and download files based on DOI.

    Parameters
    ----------
    permanent_doi : str
        Permanent DOI for a Zenodo record.
    file_name : str
        Name of the file to be downloaded.
    """

    def __init__(self, permanent_doi, file_name):
        self.permanent_doi = permanent_doi
        self.file_name = file_name
        self.headers = {
            "X-Requested-With": "XMLHttpRequest",
            "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36",
            "Accept": "*/*",
            "Host": "zenodo.org",
            "Connection": "keep-alive",
            "Pragma": "no-cache",
            "Cache-Control": "no-cache, no-store, must-revalidate",
            "Expires": "0",
        }
        self.latest_doi = self.get_latest_zenodo_doi()

    def get_latest_zenodo_doi(self):
        """Retrieves the latest version DOI using the Zenodo API."""
        url = f"https://zenodo.org/api/records/{self.permanent_doi}"
        r = requests.get(url, headers=self.headers)

        if r.status_code != 200:
            raise ValueError(f"Could not retrieve metadata: {r.status_code}")

        record_data = r.json()

        # Check if "doi" exists in the response
        latest_doi = record_data.get("doi")
        if not latest_doi:
            return self.permanent_doi  # Default to original DOI if not found

        latest_doi_number = latest_doi.split("10.5281/zenodo.")[-1]
        print(f"Latest DOI N° : {latest_doi_number}")
        return latest_doi_number

    def download_zenodo_data(self):
        """Downloads the file from Zenodo based on the provided DOI and file name."""
        try:
            r = requests.get(
                f"https://zenodo.org/api/records/{self.latest_doi}",
                headers=self.headers,
            )
            r.raise_for_status()

            record_data = r.json()
            files = record_data["files"]

            file_to_download = next(
                (f for f in files if f["key"] == self.file_name), files[0]
            )

            file_url = file_to_download["links"]["self"]

            # Check for partially downloaded file
            downloaded = 0
            if os.path.exists(self.file_name):
                downloaded = os.path.getsize(self.file_name)

            headers = self.headers.copy()
            if downloaded > 0:
                headers["Range"] = f"bytes={downloaded}-"

            headers = self.headers.copy()
            if downloaded > 0:
                headers["Range"] = f"bytes={downloaded}-"

            response = requests.get(file_url, stream=True)
            response.raise_for_status()

            file_size = int(response.headers.get("content-length", 0)) + downloaded
            mode = "ab" if downloaded > 0 else "wb"

            with open(self.file_name, mode) as file:
                with tqdm(
                    total=file_size,
                    initial=downloaded,
                    unit="B",
                    unit_scale=True,
                    desc="Downloading",
                ) as progress:
                    for chunk in response.iter_content(chunk_size=8192):
                        file.write(chunk)
                        progress.update(len(chunk))

            print(f"File downloaded successfully to {self.file_name}.")
        except requests.RequestException as e:
            print(f"Error: {e}")


zenodo_downloader = ZenodoDownloader(permanent_doi, file_name)
zenodo_downloader.download_zenodo_data()
