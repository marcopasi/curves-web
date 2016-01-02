# license
import tempfile
from urllib2 import urlopen
import os

def download_pdb(pdbid):
    """
    Download the specified PDB file in a temporary location.

    The user is responsible for deleting the file when done with it.
    """
    url='http://www.rcsb.org/pdb/files/%s.pdb'%pdbid
    response = urlopen(url)
    pdbcontent = response.read()
    pdbfile, pdbfilename = tempfile.mkstemp(suffix=".pdb") #, dir=app.static_folder)
    siz = os.write(pdbfile, pdbcontent)
    os.close(pdbfile)
    return pdbfilename
