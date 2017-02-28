# Copyright (C) 2015-2017 Marco Pasi <mf.pasi@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
"""
The Curves+ run configuration form.
"""

from flask.ext.wtf import Form
from wtforms import Form as UnsecureForm
from wtforms.fields import TextField, BooleanField, RadioField, \
     FloatField, FormField, FieldList, IntegerField, SelectField
from wtforms.validators import Required, Length, \
     ValidationError
from flask.ext.wtf.file import FileField  # , FileRequired, FileAllowed
import re


# --- validators
class Sequence(object):
    """ Validate sequence of integers, separated by ' ' or ':'. """
    sequence_re = re.compile(r"^(?:\d+(?:[: ]|$))*$")

    def __call__(self, form, field):
        if not self.sequence_re.match(field.data):
            raise ValidationError(
                'Not a valid sequence: <{}>'.format(field.data))


# --- forms
class StrandForm(UnsecureForm):
    """
    Sub-form to describe a strand: included subunits and direction.
    As it extends the non-secure Form, it can never be used alone.
    """
    nucleotides = TextField('Nucleotides', [Sequence()])
    direction = RadioField('Direction', [Required()], choices=[('+1', ''),
                                                               ('-1', '')])


class CurvesForm(Form):
    """ Form to describe a Curves+ run. """
    num_strands = 4

    # pdbfile = FileField('File', [Required()])
    # pdbid = TextField('PDBid', [Required(), Length(min=4, max=4)])
    pdbfile = FileField('File')
    pdbid = TextField('PDBid')

    strands = FieldList(FormField(StrandForm), min_entries=num_strands)

    fit = BooleanField('Fit')
    circ = BooleanField('Circ')
    line = BooleanField('Line')
    zaxe = BooleanField('Zaxe')
    refo = BooleanField('Refo')
    # test = BooleanField('Test')
    boolean_fields = 'fit circ line zaxe refo'.split()

    viewer = SelectField('3D Viewer', choices=[
        ('jsmol', 'JSmol (slow but fully featured)'),
        ('ngl', 'NGL (fast and simple)'),
        ('none', 'None (output files and plots)')])

    back = TextField('Back',  [Required(), Length(min=1, max=64)])
    wback = FloatField('Wback', [Required()])
    wbase = FloatField('Wbase', [Required()])
    naxlim = IntegerField('Naxlim', [Required()])
    rvfac = FloatField('RVFac', [Required()])

    def validate_pdbfile(self, field):
        pass
        # app.logger.info(
        #     "PDB <"+self.pdbfile.data+"> ID<"+self.pdbid.data+">")
        # if len(self.pdbfile.data) == 0 and len(self.pdbid.data) == 0:
        #     raise ValidationError('Must specify either PDB File of PDBid')

    validate_pdbid = validate_pdbfile

    def validate_strands(self, field):
        """ Accept if subunits are specified for any of the strands."""
        if not reduce(
                lambda x, y: x or y,
                [len(entry.data['nucleotides']) > 0
                 for entry in field.entries]):
            raise ValidationError(
                'Must specify subunits for at least one strand.')
