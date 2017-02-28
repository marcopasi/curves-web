# Copyright (C) 2015-2017 Marco Pasi <mf.pasi@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
from flask import Flask, request, session, redirect, \
    render_template, flash, make_response
from flask.ext.uploads import UploadSet, configure_uploads
from flask_bootstrap import Bootstrap
from flask_nav import Nav
from flask_nav.elements import Navbar, View, Subgroup

import os
from functools import wraps, update_wrapper
from datetime import datetime
from zipfile import ZipFile
import libcurves

# --- Uncomment the following lines to have lis and pdb files downloaded
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

# --- Local imports: require `app` to be defined
from .util import set_umask, download_pdb
from .forms import CurvesForm
from .curvesrun import WebCurvesConfiguration, DummyCurvesRun, \
    SubprocessCurvesRun

# --- Umask definition
# Use util.set_umask to change both os.umask and the umask
# of tempfiles created using util.mkstemp and util.mkdtemp.
# This may be useful when the web application is run
# by a user who is not www-data.
#
# Note 1: All new files respect os.umask except tempfile.*
# Note 2: shutil.copy preserves permissions
# Note 3:
#  1) the temporary input PDB file is created either:
#     - in util.py (download_pdb) using util.mkstemp
#     - here, further down (pdbfiles.save)
#  2) the temporary folder is created in curvesrun.py
#     SubprocessCurvesRun.__init__ uses util.mkdtemp
#  3) the input PDB is copied locally in curvesrun.py
#     CurvesRun.__init__ uses shutil.copy
#  4) png files are created here (plot) using pyplot.savefig
#  5) zip file is created here (makezip) using zipfile
#
set_umask(0027)

# --- Uncomment for full Debug
# from flask_debugtoolbar import DebugToolbarExtension
# toolbar = DebugToolbarExtension(app)

# --- Backend actions around request
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


# --- "nocache" wrapper
def nocache(view):
    @wraps(view)
    def no_cache(*args, **kwargs):
        response = make_response(view(*args, **kwargs))
        response.headers['Last-Modified'] = datetime.now()
        response.headers['Cache-Control'] = """\
no-store, no-cache, must-revalidate, post-check=0, pre-check=0, max-age=0"""
        response.headers['Pragma'] = 'no-cache'
        response.headers['Expires'] = '-1'
        return response
    return update_wrapper(no_cache, view)


# --- "static" routes
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


@app.route('/', methods=['GET'])
def home():
    return render_template('home.html')


# --- Success response
def success_response(run, conf=None):
    """
    Handle the response for a successful Curves+ run.
    Specify the CurvesRun object and optionally the
    CurvesConfiguration.
    """

    def OutputFile(path, url, name, ext):
        """ Utility function to generate file information dict """
        return dict(path=path,
                    url=url,
                    name=name,
                    extension=ext)

    """ Generate the =files= structure to pass to the template """
    files = {}
    # extensions = []
    for ext in run.output_extensions:
        if not os.path.isfile(run.output_file(ext)):
            files[ext] = None
        else:
            files[ext] = OutputFile(run.output_file(ext),
                                    run.output_url(ext),
                                    run.outfile+ext, ext)
    files["in"] = OutputFile(
        os.path.join(run.outdir, run.infile),
        run.urlbase+"/"+run.infile,
        run.infile, ".pdb")

    """ Fill session with necessary information from run """
    session['lisfile'] = run.output_file(".lis")
    session['outdir'] = run.outdir
    session['outurl'] = run.urlbase
    session['files'] = files

    return render_template('analyse.html',
                           files=files,
                           job=run.jobname,
                           options=conf)


# --- Analyse
@app.route('/analyse', methods=['GET', 'POST'])
def analyse():
    configuration = WebCurvesConfiguration()
    form = None
    form = CurvesForm(request.form, obj=configuration)
    # if len(request.form) == 0:
    #     form = CurvesForm(obj=configuration)
    # else:
    #     form = CurvesForm(request.form)

    class NoException(Exception):
        """Dummy exception to exit a try clause."""
        pass

    try:
        if request.method == 'GET' and request.args.get("runID") is not None:
            runid = request.args.get("runID")
            outdir = os.path.join(app.static_folder, runid)
            urlbase = outdir[outdir.find('static'):]
            curvesrun = DummyCurvesRun(outdir=outdir, urlbase=urlbase,
                                       jobname=runid)
            app.logger.info(curvesrun)
            if curvesrun.any_output():
                app.logger.info("Displaying run {}".format(runid))
                return success_response(curvesrun, configuration)
            else:
                flash("""\
No previously executed run with ID {} found.""".format(runid),
                      'danger')
                raise NoException()

        elif request.method == 'POST' and form.validate():
            "Populate configuration with form entries."
            form.populate_obj(configuration)
            # "Take some form entries as options."
            # options.update_from(form.data)
            app.logger.info("Configuration: <{}>".format(configuration))

            """Treat special case of input file"""
            jobname = None
            # XXX Taken from flask-wtf.file.FileField.has_file():
            #     is this sufficient to guarantee a valid file is available?
            if 'pdbfile' in request.files and \
               request.files['pdbfile'].filename not in [None, '', '<fdopen>']:
                app.logger.info("Saving PDB file")
                basename = pdbfiles.save(request.files['pdbfile'])
                pdbfilename = pdbfiles.path(basename)
                jobname = request.files['pdbfile'].filename
            elif len(configuration.pdbid) == 4:
                app.logger.info("Downloading PDB ID")
                pdbfilename = download_pdb(configuration.pdbid)
                jobname = configuration.pdbid
            else:
                flash(
                    "ERROR: Must specify either a PDB File or a valid PDBid.",
                    "danger")
                raise NoException()

            app.logger.info("""\
JOB name: <{}>
PDB file: <{}>""".format(jobname, pdbfilename))

            curvesrun = SubprocessCurvesRun(configuration, pdbfilename,
                                            jobname=jobname)
            retrun = curvesrun.run()      # XXX TODO: async

            message = ""
            if len(curvesrun.stderr):
                message = " Curves+ said: '{}'".format(curvesrun.stderr)

            if not curvesrun.any_output():
                """If all output files aren't present, flash stderr and stop"""
                flash("ERROR: No output files produced. "+message, 'danger')
                raise NoException()
            elif not curvesrun.all_outputs():
                """If any output file isn't present, flash stderr"""
                message = "WARNING: Some output files missing. "+message

            if len(message):
                flash(message, 'warning')
                app.logger.info("FLASH: <{}>".format(message))

            app.logger.info("Temp dir <{}>".format(curvesrun.urlbase))
            if retrun:
                return success_response(curvesrun, configuration)
    except NoException:
        pass
    except AssertionError as e:
        raise e
    except Exception as e:
        if app.config["DEBUG"]:
            import traceback
            traceback.print_exc()
        else:
            flash(
                "An error occured: Please try again. (ERROR: %s)".format(e),
                'danger')
    return render_template('prepare.html', form=form)


# --- Zip file
@app.route('/zip', methods=['GET'])
@nocache
def makezip(filename="output"):
    filename = request.args.get("prefix", filename)
    zipname = os.path.join(session['outdir'], filename+".zip")
    with ZipFile(zipname, 'w') as outzip:
        for output in session['files'].itervalues():
            outzip.write(output['path'], filename+output['extension'])
    return redirect(session['outurl']+"/"+filename+".zip")


# --- Plots
@app.route('/plot/<string:variable>', methods=['GET'])
@nocache
def plot(variable):
    if variable is not None and not libcurves.Curves.is_variable(variable):
        # flash("Variable <{}> not recognised".format(variable), 'danger')
        return "Variable <{}> not recognised".format(variable)

    filename = variable+".png"
    filepath = os.path.join(session['outdir'], filename)
    fileurl = "/"+session['outurl']+"/"+filename
    if not os.path.isfile(filepath):
        app.logger.info("Making plot <{}> <{}>".format(
            variable, datetime.now()))
        try:
            curves = libcurves.Curves(session['lisfile'])
            plot = curves.plot(variable)
            fig = plot[0].figure
            fig.gca().grid()
            # app.logger.info(filepath+":"+fileurl)
            fig.savefig(filepath)
        except:
            # flash("ERROR: Couldn't produce plot.", 'danger')
            # return analyse()
            if app.config["DEBUG"]:
                import traceback
                traceback.print_exc()
                return "ERROR: Couldn't produce plot."
    app.logger.info("Done plot <{}> <{}>".format(
        variable, datetime.now()))
    return redirect(fileurl)

# --- Navbar
nav = Nav()


@nav.navigation()
def navbar():
    return Navbar(
        View('Curves+', 'home'),
        View('Web server', 'analyse'),
        Subgroup(
            'Documentation',
            View('User Instructions', 'instructions'),
            View('Helical parameter guide', 'helpar'),
            View('Backbone parameter guide', 'bbpar')),
        View('Citing Curves+', 'cite'),
        View('Downloads', 'misc'))

nav.init_app(app)

# --- Main
if __name__ == '__main__':
    app.run(host='0.0.0.0')
