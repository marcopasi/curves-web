# Copyright (C) 2015-2016 Marco Pasi <mf.pasi@gmail.com> 
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
"""
Utilities for the Curves+ website.
"""
import tempfile
from urllib2 import urlopen
import os


class umasktempfile(object):
    """
    Define wrappers for tempfile.mkdtemp and tempfile.mkstemp
    that respect a globally set umask.
    """
    _umask = 0022

    def __init__(self, umask=None):
        if umask is not None:
            self.umask = umask

    # --- Properties
    @property
    def umask(self):
        return self._umask

    @umask.setter
    def umask(self, umask):
        self._umask = umask
        os.umask(self.umask)

    def set_umask(self, umask):
        self.umask = umask

    # --- Utilities
    def _get_permissions(self, base_):
        return base_ & (~self.umask)

    def _get_file_permissions(self):
        return self._get_permissions(0666)

    def _get_dir_permissions(self):
        return self._get_permissions(0777)

    def mkstemp(self, *args, **kwargs):
        """
        Wrapper for tempfile.mkstemp that respects the globally set umask.
        """
        file_object, file_name = tempfile.mkstemp(*args, **kwargs)
        os.chmod(file_name, self._get_file_permissions())
        return file_object, file_name

    def mkdtemp(self, *args, **kwargs):
        """
        Wrapper for tempfile.mkdtemp that respects the globally set umask.
        """
        dir_name = tempfile.mkdtemp(*args, **kwargs)
        os.chmod(dir_name, self._get_dir_permissions())
        return dir_name

tempfilewrapper = umasktempfile()

set_umask = tempfilewrapper.set_umask
mkdtemp = tempfilewrapper.mkdtemp
mkstemp = tempfilewrapper.mkstemp


def download_pdb(pdbid):
    """
    Download the specified PDB file in a temporary location.

    The user is responsible for deleting the file when done with it.
    """
    url = 'http://www.rcsb.org/pdb/files/{}.pdb'.format(pdbid)
    response = urlopen(url)
    pdbcontent = response.read()
    pdbfile, pdbfilename = mkstemp(suffix=".pdb")
    os.write(pdbfile, pdbcontent)
    os.close(pdbfile)
    return pdbfilename
