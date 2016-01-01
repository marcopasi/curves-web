### XXX license
from flask.ext.wtf import Form
from wtforms import Form as UnsecureForm
from wtforms.fields import TextField, BooleanField, FileField, \
     RadioField, FloatField, FormField, FieldList
from wtforms.validators import Required, Length, AnyOf, regexp, \
     ValidationError, StopValidation
import re

from . import app

#--- validators
class Sequence(object):
    """ Validate sequence of integers, separated by ' ' or ':'. """
    sequence_re = re.compile(r"^(?:\d+(?:[: ]|$))*$")
    def __call__(self, form, field):
        if not self.sequence_re.match(field.data):
            raise ValidationError('Not a valid sequence: <%s>'%(field.data))


class EitherRequired(Required):
    """ Stop checking field if another field is not empty. """
    def __init__(self, other_field_name, *args, **kwargs):
        self.other_field_name = other_field_name
        super(EitherRequired, self).__init__(*args, **kwargs)

    def __call__(self, form, field):
        assert self.other_field_name in form, \
               'No such field %s.'%self.other_field_name
        other_field = form[self.other_field_name]
        if bool(other_field.data):
            raise StopValidation() # no error message
        super(EitherRequired, self).__call__(form, field)
        


#--- forms
class StrandForm(UnsecureForm):
    """
    Sub-form to describe a strand: included subunits and direction.

    As it extends the non-secure Form, it can never be used alone.    
    """
    nucleotides = TextField('Nucleotides', [Sequence()])
    direction   = RadioField('Direction', [Required()], choices=[('+1',''), ('-1','')])

    
class CurvesForm(Form):
    """ Form to describe a Curves+ run. """
    num_strands = 4
        
    pdbfile = FileField('File', [EitherRequired("pdbid"), regexp('pdb$')])
    pdbid = TextField('PDBid', [EitherRequired("pdbfile"), Length(min=4, max=4)])

    strands = FieldList(FormField(StrandForm), min_entries=num_strands)

    fit = BooleanField('Fit')
    circ = BooleanField('Circ')
    line = BooleanField('Line')
    zaxe = BooleanField('Zaxe')
    refo = BooleanField('Refo')
    test = BooleanField('Test')
    boolean_fields = 'fit circ line zaxe refo test'.split()

    back  = TextField('Back',  [Required(), Length(min=1, max=64)])
    wback = FloatField('Wback', [Required()])
    wbase = FloatField('Wbase', [Required()])

    # def validate_pdbfile(self, field):
    #     if len(self.pdbfile.data) == 0 and len(self.pdbid.data) == 0:
    #         raise ValidationError('Must specify either PDB File of PDBid')
    
    # validate_pdbid = validate_pdbfile

