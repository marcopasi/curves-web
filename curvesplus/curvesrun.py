# license
import re
import os
from shutil import copy
import tempfile

from . import app

class CurvesConfiguration(object):
    """ """
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
    test = False

    back  = 'P'
    wback = 2.9
    wbase = 3.5

    naxlim = 3

    _params = """fit circ line zaxe refo test 
    back wback wbase 
    naxlim""".split()
        
    def __init__(self, config={}, **kwargs):
        self.update(config)
        self.update(kwargs)

    def update(self, dict_):
        self.__dict__.update(dict_)

    @staticmethod
    def fortran_tf(boo):
        if boo: return ".t."
        return ".f."
    
    def _param_string(self, params):
        """ Generate a valid namelist string for the specified parameters. """
        for param_name in params:
            value = getattr(self, param_name)
            if type(value) is bool:
                value = self.fortran_tf(value)
            yield "{0}={1}".format(param_name, value)

    def _strand_string(self, strands):
        """ Generate a valid strand specification string for the specified strands. """
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
        """ Generate the complete configuration string for a Curves+ run. """
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
    """
    output_extensions = ".lis .cda _X.pdb _B.pdb".split()
    
    def __init__(self, config, lib, exefile, infile, outdir, outfile="output", infile_local_name="input.pdb"):
        self.config = config
        libfile = lib+"_b.lib"
        assert os.path.isfile(libfile), "Library file doesn't exist <{}>".format(libfile)
        assert os.path.isfile(exefile), "Executable file doesn't exist <{}>".format(exefile)
        assert os.path.isfile(infile), "Input file doesn't exist <{}>".format(infile)
        self.outdir = outdir
        self.outfile = outfile
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
        return reduce(lambda x,y: x or y, 
            [os.path.isfile(output) for output in self.output_files()])

    def any_output(self):
        return self._check_outputs()
    
    def all_outputs(self):
        return self._check_outputs(op = lambda x,y: x and y)
    
    def run(self):
        pass

class SubprocessCurvesRun(CurvesRun):
    """
    Uses subprocess
    Outdir is a tempfile in the static/
    libfile and executable taken from configuration
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

        cwd = os.getcwd()
        os.chdir(self.outdir)
        p = Popen([ self.executable ],
              stdout=PIPE, stdin=PIPE, stderr=PIPE)
        [stdout, stderr] = p.communicate(self.config_string)

        self.stdout = stdout
        self.stderr = stderr
        self.returncode = p.returncode

        os.chdir(cwd)
        if p.returncode == 0:
            return True
        return False
