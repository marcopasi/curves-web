# Copyright (C) 2015-2016 Marco Pasi <mf.pasi@gmail.com> 
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
from flask import Flask, request, session, g, redirect, url_for, \
     abort, render_template, flash
from flask.ext.uploads import UploadSet, configure_uploads
from flask_bootstrap import Bootstrap, StaticCDN
from flask_nav import Nav
from flask_nav.elements import Navbar, View, Subgroup

import os
from zipfile import ZipFile
import libcurves
# import mimetypes
# mimetypes.add_type('text/plain', '.lis')
# mimetypes.add_type('text/plain', '.pdb')

app = Flask(__name__, instance_relative_config=True)
app.config.from_object('config')
app.config.from_pyfile('config.py')
Bootstrap(app)
pdbfiles = UploadSet('pdbfiles', ('pdb',))
app.config['UPLOADED_PDBFILES_DEST'] = app.static_folder
configure_uploads(app, (pdbfiles,))
#from flask_debugtoolbar import DebugToolbarExtension
#toolbar = DebugToolbarExtension(app)

from .util import download_pdb
from .forms import CurvesForm
from .curvesrun import CurvesConfiguration, SubprocessCurvesRun

## Backend actions around request
# @app.before_request
# def before_request():
#     """ Executed before request. """
#     pass

# @app.after_request
# def after_request(response):
#     """ Executed after request. """
#     return response

# @app.teardown_request
# def teardown_request(exception):
#     """ Executed after request, even in case of exception. """
#     pass


#----- "static" routes
@app.route('/instructions', methods=['GET'])
def instructions():
    return render_template('instructions.html')

@app.route('/helpar', methods=['GET'])
def helpar():
    return render_template('helpar.html')

@app.route('/bbpar', methods=['GET'])
def bbpar():
    return render_template('bbpar.html')

@app.route('/cite', methods=['GET'])
def cite():
    return render_template('cite.html')

@app.route('/misc', methods=['GET'])
def misc():
    return render_template('misc.html')


#-----
@app.route('/', methods=['GET','POST'])
def analyse():
    configuration = CurvesConfiguration()
    form = None
    if len(request.form) == 0:
        form = CurvesForm(obj=configuration)
    else:
        form = CurvesForm(request.form)

    try:
        if request.method == 'POST' and form.validate():
            form.populate_obj(configuration)
            jobname = None
            #XXX Taken from flask-wtf.file.FileField.has_file():
            #    is this sufficient to guarantee a valid file is available?
            if 'pdbfile' in request.files and request.files['pdbfile'].filename not in [None, '', '<fdopen>']:
                app.logger.info("Saving PDB file")
                basename = pdbfiles.save(request.files['pdbfile'])
                pdbfilename = pdbfiles.path(basename)
                jobname = request.files['pdbfile'].filename
            elif len(configuration.pdbid) == 4:
                app.logger.info("Downloading PDB ID")
                pdbfilename = download_pdb(configuration.pdbid)
                jobname = configuration.pdbid
            else:
                flash('ERROR: Must specify either a PDB File or a valid PDBid.', 'danger')
                raise ValueError()

            app.logger.info("JOB name: <%s>"%jobname)
            app.logger.info("PDB file: <%s>"%pdbfilename)

            curvesrun = SubprocessCurvesRun(configuration, pdbfilename)
            retrun = curvesrun.run()      #XXX TODO: async

            message=""
            if len(curvesrun.stderr):
                message =" Curves+ said: '%s'"%(curvesrun.stderr)

            if not curvesrun.any_output():
                """ if all output files aren't present, flash stderr and stop """
                flash("ERROR: No output files produced. "+message, 'danger')
                raise ValueError()
            elif not curvesrun.all_outputs():
                """ if any output file isn't present, flash stderr """
                message = "WARNING: Some output files missing. "+message

            if len(message):
                flash(message, 'warning')
                app.logger.info("FLASH: <%s>"%message)

            app.logger.info("Temp dir <%s>"%curvesrun.urlbase)
            if retrun:
                session['lisfile'] = curvesrun.output_file(".lis")
                session['outdir']  = curvesrun.outdir
                session['outurl']  = curvesrun.urlbase
                files = {}
                extensions = []
                for ext in curvesrun.output_extensions:
                    if not os.path.isfile(curvesrun.output_file(ext)): continue
                    files[ext] = OutputFile(curvesrun.output_file(ext),
                                            curvesrun.output_url(ext),
                                            curvesrun.outfile+ext, ext)
                files["in"] = OutputFile(
                    os.path.join(curvesrun.outdir, curvesrun.infile),
                    curvesrun.urlbase+"/"+curvesrun.infile,
                    curvesrun.infile, ".pdb")
                session['files'] = files
                return render_template('analyse.html', files=files, job=jobname)
    except ValueError:
        pass
    except Exception as e:
        flash("An error occured: Please try again. (ERROR: %s)"%e, 'danger')
    return render_template('prepare.html', form=form)

def OutputFile(path, url, name, ext):
    return dict(path = path,
                url = url,
                name = name,
                extension = ext)


#-----
@app.route('/zip', methods=['GET'])
def makezip(filename="output"):
    filename = request.args.get("prefix", filename)
    zipname = os.path.join(session['outdir'],filename+".zip")
    with ZipFile(zipname, 'w') as outzip:
        for output in session['files'].itervalues():
            outzip.write(output['path'], filename+output['extension'])
    return redirect(session['outurl']+"/"+filename+".zip")


#-----
@app.route('/plot/<string:variable>', methods=['GET'])
def plot(variable):
    if variable is not None and not libcurves.Curves.is_variable(variable):
        #flash("Variable <%s> not recognised"%variable, 'danger')
        return "Variable <%s> not recognised"%variable

    try:
        curves = libcurves.Curves(session['lisfile'])
        plot = curves.plot(variable)
        fig = plot[0].figure
        fig.gca().grid()
        filename = variable+".png"
        filepath = os.path.join(session['outdir'], filename)
        fileurl  = "/"+session['outurl']+"/"+filename
        #app.logger.info(filepath+":"+fileurl)
        fig.savefig(filepath)
        return redirect(fileurl)
    except:
        #flash("ERROR: Couldn't produce plot.", 'danger')
        #return analyse()
        if app.config["DEBUG"]:
            import traceback
            traceback.print_exc()
        return "ERROR: Couldn't produce plot."


#-----
nav = Nav()

@nav.navigation()
def navbar():
    return Navbar(
        'Curves+',
        View('Web server', 'analyse'),
        Subgroup(
            'Documentation',        
            View('User Instructions', 'instructions'),
            View('Helical parameter guide', 'helpar'),
            View('Backbone parameter guide', 'bbpar')),
        View('Citing Curves+', 'cite'),
        View('Downloads', 'misc'))

nav.init_app(app)
        
#-----
if __name__ == '__main__':
    app.run(host='0.0.0.0')
