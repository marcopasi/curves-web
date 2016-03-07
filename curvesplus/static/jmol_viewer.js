/*
 * Copyright (C) 2015-2016 Marco Pasi <mf.pasi@gmail.com> 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 */

/**** Jmol for Curves+
  * version 1.4
  ****/

/* GLOBAL SETTINGS */
var JMOLPATH = "/static/vendor"; // path to the Jmol installation folder
var LOADERIMG="/static/img/curves-loader.gif"; // path to ajax loader image

var APPWIDTH= "95%";	  // !! modify .jmol width in style definition
var APPHEIGHT= 400;

var DEBUG = false;

/* PDB parsing SETTINGS */
var BBRES = "BAC";
var GRRES = "GRV";
//var AXRES = "AXI";

/* frame indexes */
var axf = 1;
var bbf = 2;
var pdf = 3; // must be last

/* VISUALIZATION SETTINGS */
// PDB: nucleic acid
var nwf = 0.15; // default wireframe radius for pdb
var ncwf= 0.15; // default wireframe radius for pdb in CPK
var nlc  = "gray"; // color for lines representation
var ncpk= 0.3; // default CPK radius for pdb in CPK
var nclc = "CPK"; // color for CPK representation
var nisoc = "yellow translucent";
// PDB: protein
var pwf = 0.15; // default wireframe radius for pdb
var plc  = "steelblue"; // color for lines representation
var pcc = "deepskyblue"; // color for CPK representation
var pisoc = "deepskyblue translucent";
// Axis and backbone
var awf = 0.3; // default wireframe radius for axis
var bwf = 0.3; // default wireframe radius for backbone
// Display modes definition: nucleic acid
var PDISPMODES = [
    "isosurface piso off;wireframe off;cartoon;color "+pcc+";",
    "isosurface piso off;wireframe -"+pwf+";color "+plc+";",
    "wireframe off;cpk off;cartoon off; \
if(piso){isosurface piso on;} else {\
print 'isocreate';piso=on;\
frame last;isosurface piso resolution 1 solvent 1.4;color isosurface "+pisoc+";\
frame all;print 'isodone';};"];
var PDISPNAMES = ["Cartoon", "Lines", "Surface"];
// Display modes definition: protein
var NDISPMODES = [
    "isosurface niso off;wireframe -"+nwf+";color "+nlc+";",
    "isosurface niso off;wireframe -"+ncwf+"; cpk "+ncpk+";color "+nclc+";",
    "wireframe off;cpk off; \
if(niso){isosurface niso on;} else {\
print 'isocreate';niso=on;\
frame last;isosurface niso resolution 1 solvent 1.4;color isosurface "+nisoc+";\
frame all;print 'isodone';};"];
var NDISPNAMES = ["Lines", "CPK", "Surface"];
// colors
var GRVCOLORS=["lawngreen", "orange", "pink", "silver"];
var BACCOLOR = "red";
var BGCOLORS =["[0xd4d4d4]","white","black","powderblue"];
var AXCOLOR = "blue";

var hasBB = true;
var hasAX = true;
var hasProtein = false;
var naon = true;		// whether the na/pr are shown
var pron = true;
var naid = "NA";		// id of na/pr in Radiobutton html id's
var prid = "PR";


var loadstr = "";
/*
 * Initialize Jmol, create applet and start creating interface.
 * Also build the load script which will be run on ready.
 */
function jmol_viewer(AXPATH, BBPATH, PDBPATH) {
    jmolInitialize(JMOLPATH);
    jmolSetAppletColor(BGCOLORS[0]);
    jmolSetCallback("messageCallback", "_msg")
    jmolSetCallback("readyFunction", jmol_init)

    var loadax = "";
    var loadbb = "";
    var append = "";
    if(DEBUG) {
        loadstr += "\
console;set debug ON;\
";
    }
    /*
     * Check that the AX file is present; otherwise push the index for PDB and BB one up
     */
    if(typeof(AXPATH) != "undefined") {
        // calculate AXis rotation only if present
        loadax = "load \""+AXPATH+"\"; \
select */"+axf+"; \
cpk off; \
wireframe "+awf+"; \
color "+AXCOLOR+"; \
As={(*/"+axf+")[1]}.xyz; \
Af={(*/"+axf+")[0]}.xyz; \
df=As-Af; \
vf= @{df/sqrt(df*df)}; \
AyQ = !quaternion(cross({0 1 0},vf), acos({0 1 0}*vf));"
        append = "append";
    }else{
        hasAX = false;
        pdf -= 1;
        bbf -= 1;
    }

    /*
     * Check that the BB file is present; otherwise push the index for the PDB one up
     */
    if(typeof(BBPATH) != "undefined") {
        loadbb = "load "+append+" \""+BBPATH+"\"; \
select */"+bbf+"; \
cpk off; \
wireframe "+bwf+"; \
color gray;";
        append = "append";
    }else{
        hasBB = false;
        pdf -= 1;
    }

    /*
     * TODO: used to have automatic zoom factor;
     * but it doesn't support 100% width: solve this.
     *
    var zoomfactor = (APPHEIGHT/APPWIDTH*100); */
    var zoomfactor = 100;
    zoomstr = "\
ZF="+zoomfactor+"; \
zoom @ZF;";
    
    /*
     * PDB must be loaded last
     */
    loadstr += loadax + loadbb + "\
load "+append+" \""+PDBPATH+"\"; \
center; \
select */"+pdf+"; \
cpk off; \
wireframe off; \
delete solvent; \
sps=off; \
piso=off; \
niso=off; \
frame all; "+zoomstr+" \
print \"loaddone\";";

    /*
     * create all GUI elements
     */
    // Jmol.Info['coverImage'] = COVERIMG;
    // Jmol.Info['coverTitle'] = "Loading...";
    // Jmol.Info['deferUncover'] = true;
    jmolApplet([APPWIDTH, APPHEIGHT]);

    // show/hide models:
    jmolHtml("<div id=\"isoload\" class=\"load\"><img class=\"loader\" src=\""+LOADERIMG+"\" />Generating...</div>");
    jmolHtml("<div id=\"load\" class=\"load\"><img class=\"loader\" src=\""+LOADERIMG+"\" />Initializing 3D viewer...</div>");
    jmolHtml("<div class=\"controls row\">");
    jmolHtml("<div class=\"lcontrols col-md-6\">");
    jmolHtml("<div id=\"nadisplay\"></div>"); // will contain disptk's for nucleic acid (_initPDB)
    if(hasAX) {
        jmolHtml("<div id=\"axdisplay\">");
        jmolCheckbox(_sh(axf,1),_sh(axf,0), "Axis", "checked");
        jmolHtml("</div>");
    }
    // placeholder for dynamically derived backbone and grooves
    jmolHtml("<div id=\"bbdisplay\" >&nbsp;</div>");
    jmolHtml("</div>");
    
    // color choice and rotation buttons
    jmolHtml("<div class=\"rcontrols col-md-6\">");
    jmolHtml("<div id=\"prodisplay\"></div>"); // will contain disptk's for protein (_initPDB)
    jmolHtml("<div class=\"colors\"> Background: ");
    var colorstr = ""
    for(ci=0; ci<BGCOLORS.length; ci++) {
        var c = BGCOLORS[ci];
        colorstr += "<a href=\"#\" class=\"dummylink\"><div id=\""+c+"\" class=\"cimg\" style=\"background-color: "+_jmol2js_color(c)+"\" onclick=\"_setbg(this)\">&nbsp;</div></a>";
    }
    jmolHtml(colorstr);
    jmolHtml("</div><div class=\"buttons\">");
    jmolButton("reset;zoom @ZF;", "Reset orientation");
    jmolBr();
    if(hasAX) {
        jmolButton("reset;rotate molecular @AyQ;center;zoom @ZF;", "Align Axis");
        jmolBr();
        jmolButton("if(sps) {spin off;}else{spin @As @Af 100;};sps=!sps;", "Spin");
        jmolBr();
    }else{
        jmolButton("if(sps) {spin off;}else{spin 100;};sps=!sps;", "Spin");
        jmolBr();
    }
    jmolHtml("</div></div></div>");
}

/*
 * Run the load script and create result-dependent interface elements.
 */
function jmol_init() {
    jmolScriptWait(loadstr);
    if(DEBUG) console.log("Done loading models.");
    var chi = JSON.parse(jmolGetPropertyAsJSON("chaininfo")); // chainInfo object
    _initPDB(chi.chaininfo);
    _nash(pdf, 1);
    _prsh(pdf, 1);
    if(hasBB) {    // only init BB if the backbone file is present
	_initBB(chi.chaininfo);
    }
    if(DEBUG) console.log("Hiding loader.");
    document.getElementById("load").style.display = "none";
}


/* SHOW/HIDE FUNCTIONS */
var shstr = "display displayed or selected;";
var histr = "hide hidden or selected;";

function __sh(sh) {		// return show/hide string 
    if(sh > 0) return shstr;
    else return histr;
}

function _sh(id, sh) {		// show/hide frame id
    return "select */"+id+";"+__sh(sh);
}

function _chsh(id, ch, sh) {	// show/hide chain ch of frame id
    return "select "+ch+"/"+id+";"+__sh(sh);
}

function _psh(id, sh, kw) {	// show/hide frame id AND kw
    return "select "+kw+" AND */"+id+";"+__sh(sh);
}

function _csh(sele, sh) {	// show/hide selection sele
    return "select "+sele+";"+__sh(sh);
}

//-----
function _dosh(id, sh, sel, htid, dispmodes, isoname) {	// show/hide pdb part
    var s = "";
    if(sh == 0) {		// hide all, including isosurface
	s += _psh(id, sh, sel)+"isosurface "+isoname+" off;";
    }else{			// show checked, or first
	var checked = 0;
	
	for(var i=0; i<3; i++)
	    if(document.getElementById("jmolRadioGroup"+htid+"disp_"+i).checked) {
		checked = i;
		break;
	    }
	
	document.getElementById("jmolRadioGroup"+htid+"disp_"+checked).checked = true;
	s += _psh(id, sh, sel)+dispmodes[checked];
    }
    return s;
}

//-----
function _nash(id, sh) {	// show/hide na
    var s = "";
    if(sh == 2) naon=true;
    if(!naon) return;
    if(sh == 0) naon=false;

    s = _dosh(id, sh, "nucleic", naid, NDISPMODES, "niso")
    return jmolScript(s);
}

//-----
function _prsh(id, sh) {	// show/hide pr
    if(!hasProtein) return;
    var s = "";
    if(sh == 2) pron=true;
    if(!pron) return;
    if(sh == 0) pron=false;

    s = _dosh(id, sh, "protein", prid, PDISPMODES, "piso")
    return jmolScript(s);
}

/* Multi-plex handling */

function __atominfo(ai) { // decode atomInfo string and return [resname, chain]
    var id0 = ai.indexOf("[");
    var id1 = ai.indexOf("]");
    var id2 = ai.indexOf(":");
    return [
	ai.substring(id0+1,id1),   // ai.substring(id1+1,id2),
	ai.substring(id2+1,id2+2)];
}

function _chain_res(chi,resn) { // find chains starting with residue "resn"
    var ret = [];
    for (var ci=0; ci<chi.chains.length; ci++) {
	var tmp = __atominfo(chi.chains[ci].residues[0].atomInfo1);
	if(tmp[0]==resn)
	    ret.push(tmp[1])
    }
    return ret
}

function _bac(chi) {return _chain_res(chi, BBRES);} // return backbone chains

function _grv(chi) {return _chain_res(chi, GRRES);} // return groove chains

function _sbbb(bbchi) {		// selection for the backbone
    var bbs = _bac(bbchi);
    var Nbb = bbs.length;
    var ret = "(";
    for (var ci=0; ci<Nbb; ci++) {
	ret += ":"+bbs[ci];
	if(ci<Nbb-1) ret += " OR "
    }
    return ret+")";
}

function _sbbg(bbchi) {		// [selections] for the grooves
    var grs = _grv(bbchi);
    var Nbb = grs.length;
    var ret = [];
    for (var ci=0; ci<Nbb; ci++)
	ret.push(":"+grs[ci]);
    return ret;
}

/* HELPER FUNCTIONS */

function _setbg(el) {		// set background color of Applet
    jmolScript("background "+el.id);
    return false;
}

function _jmol2js_color(c) {
    /*
     * jmol color in rgb is: [0xAABBCC]
     * js color in rgb is    : #AABBCC
     */
    c = c.replace(/^\s+|\s+$/g, '');
    if(c[0] == "[")		//ok, it's a RGB jmol color
	c = "#"+c.slice(3,-1);
    return c;
}

function _msg(el,msg) {		// messageCallback function
    msg=""+msg;
    if(DEBUG) console.log("JMol: ",msg);
    
    if(msg.indexOf("isocreate") >= 0) { // start creating isosurface
        if(DEBUG) console.log("Showing surface loader.");
	document.getElementById("isoload").style.display = "block";
    }else if(msg.indexOf("isodone") >= 0) { // done creating isosurface
        if(DEBUG) console.log("Hiding surface loader.");
	document.getElementById("isoload").style.display = "none";
    }
    return true;
}

//---------------
function _shwrap(prep, fun) {
    if(NDISPNAMES[prep] == "Surface")
        document.getElementById("isoload").style.display = "block";
    setTimeout(fun, 10);
    return true;
}
function _nashwrap(id, prep) {
    fun = function(){ _nash(id, 1); };
    return _shwrap(prep, fun);
}
function _prshwrap(id, prep) {
    fun = function(){ _prsh(id, 1); };
    return _shwrap(prep, fun);
}

function _initPDB(chi) { // initialise dynamically-derived pdb: nucleic acid [and protein]
    if(DEBUG) console.log("initPDB");
    var nacontrols = "";
    var procontrols = "";
    jmolSetDocument(null);
    if(DEBUG) console.log("Creating NA controls...");
    // Nucleic acid
    nacontrols += jmolHtml("<div class=\"disptk\"><div class=\"nacheck\">");
    nacontrols += jmolCheckbox("javascript _nash("+pdf+",2)","javascript  _nash("+pdf+",0)", "Nucleic Acid: ", "checked");
    nacontrols += jmolHtml("</div><div class=\"naradio\">");
    pradio = [];
    for(var prep=0; prep<NDISPMODES.length; prep++) {
	var checked = 0+(prep==(0+(!hasBB)));
	pradio.push(["javascript _nashwrap("+pdf+","+prep+")", NDISPNAMES[prep], checked]);
    }
    nacontrols += jmolRadioGroup(pradio, null, "jmolRadioGroup"+naid+"disp")
    nacontrols += jmolHtml("</div></div>");

    if(DEBUG) console.log("Creating PRO controls...");
    // determine whether a protein is present
    jmolScriptWait("select protein");
    var number_models = jmolEvaluate("{selected}.size");
    if(DEBUG) console.log("Protein models: ",number_models);
    if(number_models > 0) {
	hasProtein=true;
	// Protein
	procontrols += jmolHtml("<div class=\"disptk\"><div class=\"nacheck\">");
	procontrols += jmolCheckbox("javascript _prsh("+pdf+",2)","javascript  _prsh("+pdf+",0)", "Protein: ", "checked");
	procontrols += jmolHtml("</div><div class=\"naradio\">");
	pradio = [];
	for(var prep=0; prep<PDISPMODES.length; prep++) {
	    var checked = 0;
	    pradio.push(["javascript _prshwrap("+pdf+","+prep+")", PDISPNAMES[prep], checked]);
	}
	procontrols += jmolRadioGroup(pradio, null, "jmolRadioGroup"+prid+"disp")
	procontrols += jmolHtml("</div></div>");
    }else{
	hasProtein=false;
    }
    
    if(DEBUG) console.log("DONE.");
    document.getElementById("nadisplay").innerHTML = nacontrols;
    document.getElementById("prodisplay").innerHTML = procontrols;
    return true;
}

//---------------
function _initBB(chi) { // initialise dynamically-derived backbone
    /*
     * Color bb and grooves
     */
    var grooves;
    var bbsele;
    if(chi != null) {
	var bbchi;
	if(hasAX)
	    bbchi = chi.models[1];
	else
	    bbchi = chi.models[0];
	bbsele   = _sbbb(bbchi);
	grooves = _sbbg(bbchi);
    }else{ // fall back on 2 chains
	grooves = [':C',':D'];
	bbsele = '(:A OR :B)';
    }
    
    if(DEBUG) console.log("Building BB selection...");
    var bbstr = "";
    for(gi=0; gi<grooves.length; gi++) {
	bbstr += "select */"+bbf+" AND "+grooves[gi]+"; color bonds "+GRVCOLORS[gi]+";";
    }
    bbstr += "select */"+bbf+" AND "+bbsele+"; color "+BACCOLOR+";";
    if(DEBUG) console.log("DONE: ", bbstr);
    var ret1 = jmolScript(bbstr);
    
    /*
     * Insert the controls for bb and grooves
     */
    
    if(DEBUG) console.log("Creating BB controls...");
    var bbcontrols = "";
    jmolSetDocument(null);
    var bbsele2="*/"+bbf+" AND "+bbsele;
    bbcontrols += jmolCheckbox(_csh(bbsele2, 1), _csh(bbsele2,0), "Backbone", "checked");
    bbcontrols += jmolBr();
    var glen = grooves.length;
    for(gi=0; gi<glen; gi++) {
	if(gi < glen - 1) 
	    gid = (gi+1)+""+(gi+2);
	else
	    gid = (gi+1)+"1"
	bbcontrols += jmolCheckbox(_chsh(bbf, grooves[gi], 1), _chsh(bbf, grooves[gi], 0), "Groove "+gid, "checked");
	bbcontrols += jmolBr();
    }
    if(DEBUG) console.log("DONE");
    document.getElementById("bbdisplay").innerHTML = bbcontrols;
    return ret1;
}
