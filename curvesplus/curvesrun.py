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
from shutil import copy as shell_copy
from zipfile import ZipFile
import tempfile

from . import app

#---
class Configuration(object):
    """Flexible dict-like object for configuration."""
    def __init__(self, config={}, **kwargs):
        """ Initialise with parameters.
        Parameters can be specified both using the
        =config= option (dict of parameters) or directly
        in the constructor call."""
        self.update(config)
        self.update(kwargs)

    def update(self, dict_):
        """Update the class with the specified parameters."""
        self.__dict__.update(dict_)

    def update_existing(self, dict_):
        """Update the class with specified parameters only
        if they are alreay defined as attributes."""
        for key, value in dict_.iteritems():
            if hasattr(self, key):
                setattr(self, key, value)

    def __str__(self):
        """Stringify meaningful, dict-like information."""
        return ", ".join(
            [ "%s: %s"%(attr, getattr(self, attr))
               for attr in dir(self)
            if not attr.startswith('__') and not callable(getattr(self,attr))])


#---
class CurvesConfiguration(Configuration):
    """ Configuration of a generic Curves+ run. """
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
    rvfac = 7.5

    "These parameters go in the namelist"
    _params = """fit circ line zaxe refo 
    back wback wbase 
    naxlim
    rvfac""".split()

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


#---
class WebCurvesConfiguration(CurvesConfiguration):
    """ Configuration of a web-site Curves+ run. """
    viewer = True


#---
class DummyCurvesRun(object):
    """
    Not a Curves+ run. 
    """
    output_extensions = ".lis _X.pdb _B.pdb _R.pdb".split()
    
    def __init__(self, outdir="", config=None, libbase="", exefile="", infile="input.pdb", 
                 urlbase="", outfile="output", jobname="Curves+"):
        """ Initialise the curves run providing:
        =outdir=:  Full path to run output directory.
        =config=:  CurvesConfiguration containing parameters for the run,
        =libbase=: Base name for the Curves+ library,
        =exefile=: Full path to Curves+ executable,
        =infile=:  Name of run input file,

        The default base name of output files is "output", but can be
        configured using =outfile=. Optionally a jobname can be provided.
        """
        self.outdir = outdir
        self.config = config
        self.libbase = libbase
        self.executable = exefile
        self.infile = infile
        self.urlbase = urlbase
        self.outfile = outfile
        self.jobname = jobname

    def output_url(self, extension):
        """Return the URL of the output file with the specified extension"""
        return self.urlbase+"/"+self.outfile+extension

    def output_file(self, extension):
        """Return the PATH of the output file with the specified extension"""
        return os.path.join(self.outdir,self.outfile+extension)

    def output_files(self):
        """Return the PATH of all output files"""
        for ext in self.output_extensions:
            yield self.output_file(ext)

    def _check_outputs(self, op = lambda x,y: x or y):
        """ Return a dictionary specifiying if each output
        file exists. """
        return {output: os.path.isfile(output) for output in self.output_files()}

    def any_output(self):
        """ Shortcut to check if any output file exists. """
        return reduce(lambda x,y: x or y, self._check_outputs().values())

    def all_outputs(self):
        """ Shortcut to check if all output files exist. """
        return reduce(lambda x,y: x and y, self._check_outputs().values())


#---
class CurvesRun(DummyCurvesRun):
    """
    Abstract Curves+ run. Subclass this to generate
    useable classes by implementing the "run()" method.
    """
    def __init__(self, outdir="", libbase="", exefile="", infile="", 
                 infile_local_name="input.pdb", **kwargs):
        """ See DummyCurvesRun. Here, some checks are performed on the input.
        The run input file is copied to a
        local location before execution, by default named "input.pdb"
        configurable with =infile_local_name=.
        """
        
        """ Check that required input files exist """
        libfile = libbase+"_b.lib"
        assert os.path.isfile(libfile), "Library file doesn't exist <{}>".format(libfile)
        assert os.path.isfile(exefile), "Executable file doesn't exist <{}>".format(exefile)
        assert os.path.isfile(infile), "Input file doesn't exist <{}>".format(infile)
        """ Copy the input file locally, then remove original """
        infile_local = os.path.join(outdir, infile_local_name)
        assert not os.path.isfile(infile_local), "Output folder already contains an input file <{}>".format(infile_local)
        shell_copy(infile, infile_local)
        os.remove(infile)
        super(CurvesRun, self).__init__(
            outdir=outdir, libbase=libbase, exefile=exefile, infile=infile_local_name,
            **kwargs)
        """ Generate the run configuration string using CurvesConfiguration.string() """
        self.config_string = self.config.string(self.infile, self.outfile, self.libbase)

    def run(self):
        pass


#---
class SubprocessCurvesRun(CurvesRun):
    """
    Curves Run that uses Python's =subprocess=.
    The output folder (=outdir=) is a tempfile (generated
    using Python's =tempfile=) in the static/ folder of the
    website. The location of the "libfile" and of the Curves+
    executable are taken from configuration variables
    CURVESPLUS_HOME and CURVESPLUS_EXE, respectively.
    """
    def __init__(self, config, infile, **kwargs):
        libbase = os.path.join(app.config["CURVESPLUS_HOME"], "standard")
        exefile = app.config["CURVESPLUS_EXE"]
        outdir = tempfile.mkdtemp(dir=app.static_folder)
        urlbase = outdir[outdir.find('static'):]
        super(SubprocessCurvesRun, self).__init__(
            outdir=outdir, config=config, libbase=libbase,
            exefile=exefile, infile=infile, urlbase=urlbase,
            **kwargs)
        
    def run(self):
        from subprocess import Popen, PIPE, STDOUT
        # if output files exist, remove them
        for output in self.output_files():
            if os.path.isfile(output):
                os.remove(output)

        app.logger.info("Config string<%s>"%self.config_string)

        cwd = os.getcwd()
        os.chdir(self.outdir)
        p = Popen([ self.executable ],
              stdout=PIPE, stdin=PIPE, stderr=PIPE)
        [stdout, stderr] = p.communicate(self.config_string)

        self.stdout = stdout
        self.stderr = stderr
        self.returncode = p.returncode
        app.logger.info("Curves+: out<%s> err<%s> ret<%d>"%(stdout,stderr,p.returncode))

        os.chdir(cwd)
        if p.returncode == 0:
            return True
        return False
