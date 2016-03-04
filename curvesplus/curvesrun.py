# Copyright (C) 2015-2016 Marco Pasi <mf.pasi@gmail.com> 
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
"""
CurvesRun: Configure and run Curves+.
"""
import re
import os
from shutil import copy
from zipfile import ZipFile
import tempfile

from . import app

class CurvesConfiguration(object):
    """ Configuration of a Curves+ run. """
    pdbfile = ''
    pdbid = ''
    
    class Strand(object):
        def __init__(self, nucl = '', dire ='+1'):
            self.nucleotides = nucl
            self.direction = dire
            
    strands = [
        Strand('','+1'),
        Strand('','-1'),
        Strand('','-1'),
        Strand('','-1')
        ]

    fit = True
    circ = False
    line = False
    zaxe = False
    refo = False
    #test = False

    back  = 'P'
    wback = 2.9
    wbase = 3.5
    naxlim = 3

    "These parameters go in the namelist"
    _params = """fit circ line zaxe refo 
    back wback wbase 
    naxlim""".split()

    def __init__(self, config={}, **kwargs):
        """ Initialise the run with parameters.
        Parameters can be specified both using the
        =config= option (dict of parameters) or directly
        in the constructor call."""
        self.update(config)
        self.update(kwargs)

    def update(self, dict_):
        """Update the class with the specified parameters."""
        self.__dict__.update(dict_)

    @staticmethod
    def fortran_tf(boo):
        if boo: return ".t."
        return ".f."
    
    def _param_string(self, params):
        """ Generate a valid Curves+ namelist string for
        the specified parameters. The string can then be
        passed to the Curves+ executable."""
        for param_name in params:
            value = getattr(self, param_name)
            if type(value) is bool:
                value = self.fortran_tf(value)
            yield "{0}={1}".format(param_name, value)

    def _strand_string(self, strands):
        """ Generate a valid strand specification string for
        the specified strands. The string can then be passed
        to the Curves+ executable."""
        strandlines = []
        directions = []
        for strand in strands:
            if len(strand.nucleotides) > 0:
                strandlines.append(strand.nucleotides)
                directions.append(strand.direction)
            else:
                directions.append(0)

        summaryline = "{0} {1}".format(
            len(strandlines),
            " ".join(("{:d}".format(int(d)) for d in directions)))
        return [summaryline] + strandlines

    def string(self, infile, outfile, libprefix):
        """ Generate the complete configuration string for
         a Curves+ run, by combining the namelist string and
         the strand-specification string."""
        param_string = ", ".join(self._param_string(self._params))
        param_string = "\n".join(re.findall("(.{,64}(?:,|$))", param_string)[:-1])
        strand_string = "\n".join(self._strand_string(self.strands))
        return """\
&inp file={0}, lis={1},
 {2},
 lib={3}
&end
{4}
""".format(infile, outfile, param_string, libprefix,
           strand_string)

class CurvesRun(object):
    """
    Abstract Curves+ run. Subclass this to generate
    useable classes by implementing the "run()" method.
    """
    output_extensions = ".lis _X.pdb _B.pdb".split()
    
    def __init__(self, config, lib, exefile, infile, outdir, outfile="output", infile_local_name="input.pdb"):
        self.config = config
        libfile = lib+"_b.lib"
        assert os.path.isfile(libfile), "Library file doesn't exist <{}>".format(libfile)
        assert os.path.isfile(exefile), "Executable file doesn't exist <{}>".format(exefile)
        assert os.path.isfile(infile), "Input file doesn't exist <{}>".format(infile)
        self.outdir = outdir
        self.outfile = outfile
        self.infile = infile_local_name
        infile_local = os.path.join(self.outdir, infile_local_name)
        assert not os.path.isfile(infile_local), "Output folder already contains an input file <{}>".format(infile_local)
        copy(infile, infile_local)
        os.remove(infile)
        self.config_string = self.config.string(infile_local_name, outfile, lib)
        self.executable = exefile
        self.urlbase = ""

    def output_url(self, extension):
        return self.urlbase+"/"+self.outfile+extension

    def output_file(self, extension):
        return os.path.join(self.outdir,self.outfile+extension)

    def output_files(self):
        for ext in self.output_extensions:
            yield self.output_file(ext)

    def _check_outputs(self, op = lambda x,y: x or y):
        """ Check if output files exist. """
        return reduce(op, 
            [os.path.isfile(output) for output in self.output_files()])

    def any_output(self):
        return self._check_outputs()

    def all_outputs(self):
        return self._check_outputs(op = lambda x,y: x and y)

    def run(self):
        pass

class SubprocessCurvesRun(CurvesRun):
    """
    Curves Run that uses Python's =subprocess=.
    The output folder (=outdir=) is a tempfile (generated
    using Python's =tempfile=) in the static/ folder of the
    website. The location of the "libfile" and of the Curves+
    executable are taken from configuration variables
    CURVESPLUS_HOME and CURVESPLUS_EXE, respectively.
    """
    def __init__(self, config, infile):
        libbase = os.path.join(app.config["CURVESPLUS_HOME"], "standard")
        exefile = app.config["CURVESPLUS_EXE"]
        outdir = tempfile.mkdtemp(dir=app.static_folder)
        super(SubprocessCurvesRun, self).__init__(
            config, libbase, exefile, infile, outdir)
        self.urlbase = outdir[outdir.find('static'):]
        
    def run(self):
        from subprocess import Popen, PIPE, STDOUT
        # if output files exist, remove them
        for output in self.output_files():
            if os.path.isfile(output):
                os.remove(output)

        app.logger.info(self.config_string)

        cwd = os.getcwd()
        os.chdir(self.outdir)
        p = Popen([ self.executable ],
              stdout=PIPE, stdin=PIPE, stderr=PIPE)
        [stdout, stderr] = p.communicate(self.config_string)

        self.stdout = stdout
        self.stderr = stderr
        self.returncode = p.returncode
        app.logger.info("Curves+: %s %s %d"%(stdout,stderr,p.returncode))

        os.chdir(cwd)
        if p.returncode == 0:
            return True
        return False
