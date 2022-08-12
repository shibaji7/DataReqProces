#!/usr/bin/env python

"""zenodo.py: module is dedicated to upload files to Zenodo and create a dataset."""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2020, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import requests
import json
import os
import glob

class Zenodo(object):
    """
    This class interacts with zenondo using REST API.
    1. Create a project (fetch if exisis)
    2. Upload dataset into project
    3. Create a DOI and share
    """
    
    def __init__(self, param):
        """
        param: Name of the parameter file
        """
        self.param = param
        self.json = f"params/{param}.json"
        self.dataset = f"tmp/{param}" 
        self.ACCESS_TOKEN = "RGVC19kv3oOTMDfL9ZTuTsK1wiqKchw5TAaMOW7G8JE7wrXwpkPpK0VPKZwn"
        self.__setup__()
        self.create_proj()
        self.upload_files()
        self.publish()
        return
    
    def __setup__(self):
        self.headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self.ACCESS_TOKEN}"
        }
        self.baseurl = f"https://zenodo.org/api"
        return

    def create_proj(self):
        """
        Create an empty project
        """
        metadata = {
            "metadata": {
                "title": f"Zenodo.SD.VT.{self.param}",
                "upload_type": "dataset",
                "publication_type": "other",
                "prereserve_doi": True,
                "description": "SuperDARN dataset",
                "creators": [{"name": "Chakraborty, Shibaji", "affiliation": "SuperDARN, VT"}]
            }
        }
        r = requests.post(
            f"{self.baseurl}/deposit/depositions", 
            headers=self.headers, 
            data=json.dumps(metadata)
        )
        self.deposition_id = r.json()["id"]
        self.bucket_link = r.json()["links"]["bucket"]
        return
    
    def upload_files(self):
        """
        Upload data and files
        """
        params = {
            "access_token": self.ACCESS_TOKEN
        }
        files = glob.glob(self.dataset + "/*")
        files.append(self.json)
        for f in files:
            with open(f, "rb") as fp:
                filename = f.split("/")[-1]
                r = requests.put(
                    f"{self.bucket_link}/{filename}",
                    params=params,
                    data=fp,
                )
        return
    
    def publish(self):
        """
        Create dataset with DOI
        """
        params = {
            "access_token": self.ACCESS_TOKEN
        }
        r = requests.post(
            f"{self.baseurl}/deposit/depositions/{self.deposition_id}/actions/publish",
            params=params,
        )
        return