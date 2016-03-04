# Copyright (C) 2015-2016 Marco Pasi <mf.pasi@gmail.com> 
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
"""
LibCurves: Parse Curves+ lis outputs and plot results using matplotlib.
"""

import pandas as pd
from numpy import inf

# TODO:
# 1) FileIntervals: rename and comment
# 2) Curves: comment, cleanup, add (robuts) parsing of input parameters
# 3) LisFile: add options to read header from section
# 4) Curves: consider having explicitly 3+ series, interrelated:
#    1. Levels   (integer index of level, first column in all stanzas)
#    2. Seq  # (1 letter per level), Level direction (second, fourth column of bp_axis)
#    3. Resn # (1 number per level), Level direction (third, fifth column of bp_axis)
#    Where # indicates there is one series per strand.
# 5) LisFile: Define more robust parsing of nml entries.

#----------------------------------------------------------------------------------
class FileIntervals(file):
    """
    Raise stopIteration every time a condition is fulfilled.
    Conditions are expressed in the form of strings, which 
    must match at the beginning of the line (str.startswith).
    First section is ignored, last section is returned. To
    retain the first section, initialise with an empty string
    as the first condition.

    The sectioning feature of this class will work exclusively
    when accessing data through iterator methods. The "read"
    methods will work as usual (see file).
    """
    def __init__(self, conditions, start_deltas = [], *args, **kwargs):
        super(FileIntervals, self).__init__(*args, **kwargs)
        self.conditions = conditions
        if len(start_deltas) == 0:
            self.start_deltas = [0 for condition in conditions]
        else:
            self.start_deltas = start_deltas

        self._nc = len(self.conditions)
        assert len(self.start_deltas) == self._nc
        assert [x >= 0 for x in self.start_deltas]
        self._valid = -inf
        self._stack = []
        # initialise to first condition
        self._ic = 0

    def match(self, line, index):
        "Return true if the line matches the condition at index."
        return line.startswith(self.conditions[index])

    def next(self):
        while True:
            if len(self._stack) > 0:
                return self._stack.pop(0)
            line = super(FileIntervals, self).next()
            self._valid += 1
            if self._ic < self._nc and self.match(line, self._ic):
                delta = self.start_deltas[self._ic]
                self._ic += 1
                exvalid, self._valid = self._valid, -delta
                if exvalid >= 0:
                    if delta == 0:
                        self._stack.append(line)
                    raise StopIteration
            if self._valid >= 0:
                return line

#----------------------------------------------------------------------------------
class LisFile(FileIntervals):
    """Generic Curves/Canal LisFile class."""
    
    def __init__(self, *args, **kwargs):
        super(LisFile, self).__init__(*args, **kwargs)
    
    def read_section(self, widths, usecols, names, header=0, footer=0):
        """ Read a section of the lis file with uniform fixed-width formatting
        in a pd.DataFrame. """
        return pd.read_fwf(self, header=None, skiprows=header, skipfooter=footer,
                            widths = widths, usecols = usecols, names = names,
                            na_values=["---", "----"])

    def read_nml(self, line, keylen=7, valuelen=32, initspace=2):
        """
        Read nml lines of a lis file.

        Default line type is char, as defined in Curves+ (v2.8).
        Other types are defined as (in Curves+ v2.8):
        
        | Type  | keylen | valuelen |
        |-------+--------+----------|
        | char  | 7      | 32       |
        | int   | 7      | 6        |
        | float | 7      | 6        |
        | bool  | 7      | 6        |
        
        N.B.: Keylen includes ":".
        """
        ret = {}

        def _consume(s,a,b):
            return s[b:], s[a:b]
        
        while len(line) > 0:
            line, key  = _consume(line, initspace, initspace+keylen)
            line, value= _consume(line, 0, valuelen)
            key = key.strip()
            if ret.has_key(key):
                raise ValueError("NML Error: Duplicate configuration entry %s."%key)
            ret[key] = value.strip()
        return ret

#----------------------------------------------------------------------------------
class Curves(object):
    """
    """
    section_headers = [
        "",
        "  (A) BP-Axis", "  (B) Intra-BP parameters", 
        "  (C) Inter-BP", "  (D) Backbone Parameters",
        "  (E) Groove parameters", "  (I) " ]
    bp_axis_variables = "xdisp ydisp inclin tip ax-bend".split()
    intra_bp_variables= "shear stretch stagger buckle propeller opening".split()
    inter_bp_variables= "shift slide rise tilt roll twist h-rise h-twist".split()
    backbone_variables= "alpha beta gamma delta epsilon zeta chi phase amplitude pucker".split()
    maxgrooves = 4
    groove_variables = ["%s%d"%(var,num) for num in range(maxgrooves) for var in "width depth".split()]
    curves_variables  = bp_axis_variables + intra_bp_variables + \
                          inter_bp_variables + backbone_variables + groove_variables
    units = ("angstrom angstrom degree degree degree " +
             "angstrom angstrom angstrom degree degree degree " +
             "angstrom angstrom angstrom degree degree degree angstrom degree " +
             "degree degree degree degree degree degree degree degree degree degree").split() + \
            "angstrom angstrom".split()*maxgrooves
    assert len(units) == len(curves_variables)
    
    @classmethod
    def is_variable(cls, var):
        return var in cls.curves_variables
                              
    @classmethod
    def get_unit(cls, var):
        if not cls.is_variable(var):
            return None
        return cls.units[cls.curves_variables.index(var)]
    
    def __init__(self, filename):
        self.lisfile = LisFile(name = filename, conditions = self.section_headers)
        self._parse_lis()
        self.lisfile.close()

    def _parse_lis(self):
        """
        Read all sections of the lis file, using =section_headers=
        to define sections.
        """
        self.header = self._read_header()
        self.bp_axis = self._read_bp_axis()
        self.nbp = len(self.bp_axis)
        self.intra_bp = self._read_intrabp()
        self.inter_bp = self._read_interbp()
        self.backbone = self._read_bckbone()
        self.groove   = self._read_groove()

    def _read_header(self):
        """
        Read the Lis file header to get run parameters.
        
        These information are contained between two pairs of
        newlines.
        """
        newlines = 0
        parsing  = False
        for line in self.lisfile: # first section
            if len(line.strip()) == 0: # newline
                newlines += 1
            if newlines == 2:
                if parsing:
                    break
                parsing = True
        
    def _read_bp_axis(self): 
        #   (A) BP-Axis        Xdisp   Ydisp   Inclin    Tip  Ax-bend
        #
        #       write(6,30) i,na(i,1),nu(i,1),na(i,2),nu(i,2),var(1),var(2),
        #      1    var(4),var(5),axang
        # 30        format(2x,i3,') ',a1,i4,'-',a1,i4,2f8.2,2f8.1,a8)
        widths = (2,3,2,1,4,1,1,4) + (8,)*5
        usecols= [1,3,4,6,7] + range(8,13)
        names  = "Entry W Wn C Cn".split() + self.bp_axis_variables
        return self.lisfile.read_section(widths, usecols, names, header=2, footer=3)
        
    def _read_intrabp(self):
        #   Strands 1-2       Shear  Stretch Stagger  Buckle  Propel Opening
        #
        #       write(6,10) i,na(i,1),nu(i,1),na(i,2),nu(i,2),(var(j),j=1,6)
        # 10    format(2x,i3,') ',a1,i4,'-',a1,i4,3f8.2,3f8.1,2x,f8.2,f8.1)
        widths = (2,3,2,1,4,1,1,4) + (8,)*6           # extra fields unkown
        usecols= [1,3,4,6,7] + range(8,14)
        names  = "Entry W Wn C Cn".split() + self.intra_bp_variables
        return self.lisfile.read_section(widths, usecols, names, header=4, footer=3)
    
    def _read_interbp(self):
        #   (C) Inter-BP       Shift   Slide    Rise    Tilt    Roll   Twist   H-Ris   H-Twi
        #
        #       write(6,20) i-1,na(il,1),nu(il,1),na(iu,1),nu(iu,1),
        #      1 (var(j),j=1,6),vaf(il,1),vaf(il,2)
        # 20    format(2x,i3,') ',a1,i4,'/',a1,i4,3f8.2,3f8.1,f8.2,f8.1)
        widths = (2,3,2,1,4,1,1,4) + (8,)*8
        usecols= [1,3,4,6,7] + range(8,16)
        names  = "Entry L Ln U Un".split() + self.inter_bp_variables
        return self.lisfile.read_section(widths, usecols, names, header=2, footer=3)

    def _read_bckbone(self):
        #    Strand 1     Alpha  Beta   Gamma  Delta  Epsil  Zeta   Chi    Phase  Ampli  Puckr
        #
        #       write(6,26) i,na(i,k),nu(i,k),(tor(m),m=ms,me)
        # 26    format(2x,i3,') ',a1,i4,2x,10a7)
        widths = (2,3,2,1,4,2) + (7,)*10
        usecols= [1,3,4] + range(6,16)
        names  = "Entry L Ln".split() + self.backbone_variables
        both = self.lisfile.read_section(widths, usecols, names, header=4, footer=1)
        return [both[:self.nbp], both[self.nbp+3:]]
        
    def _read_groove(self):
        #   Level           W12     D12     W21     D21
        #
        #       write(6,40) r,string,(stg(kop,m),m=1,2*nst)
        # 40    format(2x,f5.1,a7,8(2x,a6))
        widths = (2,5,2,1,4) + (8,)*8
        usecols= [1,3,4] + range(5,13)
        names  = "N L Ln".split() + self.groove_variables
        return self.lisfile.read_section(widths, usecols, names, header=4, footer=0)

    def normalize_interval(self, interval=None):
        """Return a valid interval for this DNA molecule."""
        if interval is None:
            interval = [0, self.nbp]
        start, end = interval[:]
        # if None in interval, substitute with default
        if start is None or start < 0:
            start = 0
        if end is None or end > self.nbp:
            end = 0
        # interval[1] can be <0
        if end <= 0:
            end += self.nbp

        return start, end
        
    def get_variable(self, varname, interval=None, step=1):
        """Return the Series corresponding to an interval of the specified variable"""
        if not self.is_variable(varname):
            raise ValueError("Variable '%s' not recognised."%varname)

        start, end = self.normalize_interval(interval)
        
        if varname in self.bp_axis_variables:
            y = self.bp_axis[varname]
            x = self.bp_axis.Wn
        elif varname in self.intra_bp_variables:
            y = self.intra_bp[varname]
            x = self.intra_bp.Wn
        elif varname in self.inter_bp_variables:
            y = self.inter_bp[varname]
            x = self.inter_bp.Ln
        elif varname in self.backbone_variables:
            y = [bb[varname] for bb in self.backbone]
            x = [self.backbone.Ln, self.backbone.Ln]
        elif varname in self.groove_variables:
            y = self.groove[varname]
            x = self.groove.N
            # grooves start at 1.5, step 0.5
            start = max((start-1)*2+1, 0)
            end = (end-1)*2
        return y[start:end:step], x[start:end:step]

    @property
    def sequence(self):
        return pd.Series(self.bp_axis.W.values, index=self.bp_axis.Wn)
    @sequence.setter
    def sequence(self, sequence):
        raise NotImplementedError

    def plot(self, varname, interval=None, step=1, hold=False):
        """Plot a variable on an interval using pyplot."""
        import matplotlib.pyplot as plt

        # if start is not None and start < 1:
        #     raise ValueError("Interval start <%d> not valid."%start)
        # if end is not None and end > dna.nbp:
        #     raise ValueError("Interval end <%d> not valid (%d bp)."%(end, self.nbp))
        start, end = self.normalize_interval(interval)
        
        y, x   = self.get_variable(varname, interval = interval)
        if not hold:
            plt.clf()
        # plot inter and backbone variables at junctions
        dx = 0.0
        if varname in self.inter_bp_variables or varname in self.backbone_variables:
            dx = 0.5
        plot = plt.plot(x.values+dx, y.values)
        if dx > 0:
            # TODO: instead of extrapolating, get next value from bp_axis or intra_bp
            nx = len(x)
            nextx = 2*x.iloc[nx-1]-x.iloc[nx-2]
            if nextx in self.bp_axis.Wn.values:
                x.set_value(nx, nextx)
        # consider that grooves have 0.5 step
        if varname in self.groove_variables:
            I = x.apply(lambda x: x == round(x))
            x = x[I]
            xnames = self.groove.Ln[I].values
        else:
            xnames = x.values
        seq = self.sequence[xnames]
        if len(x) <= 20:
            plt.xticks(x, seq)
        plt.xlim([min(x), max(x)])
        if varname in self.groove_variables:
            plt.xlabel("Level")
        else:
            plt.xlabel("Basepair")
        plt.ylabel("%s (%ss)"%(varname.capitalize(), self.get_unit(varname)))
        return plot
                    
def main():
    c=Curves("/tmp/CCBD.out.lis")
    print c.backbone[1]

# ------------------------------------------
if __name__ == "__main__": main()
