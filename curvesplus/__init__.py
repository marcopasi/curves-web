# Copyright (C) 2015-2016 Marco Pasi <mf.pasi@gmail.com> 
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash
from flask.ext.uploads import UploadSet
from flask_debugtoolbar import DebugToolbarExtension
from flask_bootstrap import Bootstrap

import os
import libcurves

# configuration: see config.py
# create app
app = Flask(__name__, instance_relative_config=True)
app.config.from_object('config')
app.config.from_pyfile('config.py')
Bootstrap(app)
## app.config.from_envvar('FLASKR_SETTINGS', silent=True)
pdbs = UploadSet('pdb')
toolbar = DebugToolbarExtension(app)

from .util import download_pdb
from .forms import CurvesForm
from .curvesrun import CurvesConfiguration, SubprocessCurvesRun

# Backend actions around request
@app.before_request
def before_request():
    """ Executed before request. """
    pass

@app.after_request
def after_request(response):
    """ Executed after request. """
    return response

@app.teardown_request
def teardown_request(exception):
    """ Executed after request, even in case of exception. """
    pass

#@app.route('/analyse', methods=['GET','POST'])
@app.route('/', methods=['GET','POST'])
def analyse():
    #TODO: plots
    #TODO: 3D
    app.logger.info(request.files)
    configuration = CurvesConfiguration()
    form = None
    if len(request.form) == 0:
        form = CurvesForm(obj=configuration)
    else:
        form = CurvesForm(request.form)
        q()
        
    if request.method == 'POST' and form.validate():
        form.populate_obj(configuration)
        if 'pdbfile' in request.files:
            app.logger.info("Found PDB file")
            pdbfilename = pdbs.save(request.files['pdbfile'])
        else:
            app.logger.info("Downloading PDB ID")
            pdbfilename = download_pdb(configuration.pdbid)

        curvesrun = SubprocessCurvesRun(configuration, pdbfilename)
        retrun = curvesrun.run()      #XXX TODO: async

        message=""
        # if all output files aren't present, flash stderr
        if not curvesrun.any_output():
            message = "ERROR: No output files produced."
            retrun = False
        # if any output file isn't present, flash stderr
        elif not curvesrun.all_outputs():
            message = "WARNING: Some output files missing."
        elif len(curvesrun.stderr):
            message = "WARNING:"
            
        if len(curvesrun.stderr):
            message +=" Curves+ said: '%s'"%(curvesrun.stderr)
            
        if len(message):
            flash(message)
            app.logger.info("FLASH: <%s>"%message)
        
        if retrun:
            session['lisfile'] = curvesrun.output_file(".lis")
            session['outdir']  = curvesrun.outdir
            session['outurl']  = curvesrun.urlbase
            files = [(curvesrun.output_url(ext),
                      curvesrun.outfile+ext) for ext in curvesrun.output_extensions]
            return render_template('analyse.html', files=files)
    return render_template('prepare.html', form=form)

@app.route('/plot/<string:variable>', methods=['GET'])
def plot(variable):
    message = ""
    if variable is not None and not libcurves.Curves.is_variable(variable):
        message = "Variable <%s> not recognised"%variable

    if True:
        curves = libcurves.Curves(session['lisfile'])
        plot = curves.plot(variable)
        fig = plot[0].figure
        filename = variable+".png"
        filepath = os.path.join(session['outdir'], filename)
        fileurl  = "/"+session['outurl']+"/"+filename
        app.logger.info(filepath+":"+fileurl)
        fig.savefig(filepath)
        return redirect(fileurl)
    # except ValueError as e:


#-----
if __name__ == '__main__':
    app.run(host='0.0.0.0')
