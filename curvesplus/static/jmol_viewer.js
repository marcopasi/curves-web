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
var JMOLPATH = "/static/media"; // path to the Jmol installation folder
var LOADERIMG="/static/img/curves-loader.gif"; // path to ajax loader image

var APPWIDTH= "90%";	  // !! modify .jmol width in style definition
var APPHEIGHT= 300;

/* PDB parsing SETTINGS */
var BBRES = "BAC";
var GRRES = "GRV";
//var AXRES = "AXI";

/* frame indexes */
var axf = 1;
var bbf = 2;
var pdf = 3; // XXX must be last

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
print 'isocreate';\
frame last;isosurface piso solvent 1.4;color isosurface "+pisoc+";\
frame all;piso=on;print 'isodone';};"];
var PDISPNAMES = ["Cartoon", "Lines", "Surface"];
// Display modes definition: protein
var NDISPMODES = [
    "isosurface niso off;wireframe -"+nwf+";color "+nlc+";",
    "isosurface niso off;wireframe -"+ncwf+"; cpk "+ncpk+";color "+nclc+";",
    "wireframe off;cpk off; \
if(niso){isosurface niso on;} else {\
print 'isocreate';\
frame last;isosurface niso solvent 1.4;color isosurface "+nisoc+";\
frame all;niso=on;print 'isodone';};"];
var NDISPNAMES = ["Lines", "CPK", "Surface"];
// colors
var GRVCOLORS=["lawngreen", "orange", "pink", "silver"];
var BACCOLOR = "red";
var BGCOLORS =["[0xd4d4d4]","white","black","powderblue"];
var AXCOLOR = "blue";

var hasBB = true;
var hasAX = true;
var hasProtein = false;
var append = "";
var naon = true;		// whether the na/pr are shown
var pron = true;
var naid = "NA";		// id of na/pr in Radiobutton html id's
var prid = "PR";

/*
 * Initialize Jmol, create applet and start creating interface
 */
function jmol_viewer(AXPATH, BBPATH, PDBPATH) {

    jmolInitialize(JMOLPATH);
    jmolSetAppletColor(BGCOLORS[0]);
    jmolSetCallback("messageCallback", "_msg")

    var loadstr = "";
    /*
     * Check that the AX file is present; otherwise push the index for PDB and BB one up
     */
    if(typeof(AXPATH) != "undefined") {
        loadstr += "load \""+AXPATH+"\";"
        append = "APPEND";
    }else{
        hasAX = false;
        pdf -= 1;
        bbf -= 1;
    }

    /*
     * Check that the BB file is present; otherwise push the index for the PDB one up
     */
    if(typeof(BBPATH) != "undefined") {
        loadstr += "load "+append+" \""+BBPATH+"\";";
        append = "APPEND";
    }else{
        hasBB = false;
        pdf -= 1;
    }

    /*
     * PDB must be loaded last
     */
    
    /*
     * TODO: used to have automatic zoom factor;
     * but it doesn't support 100% width: solve this.
ZF="+(APPHEIGHT/APPWIDTH*100)+"; \
zoom @ZF; \
    */
    
    loadstr += "\
load "+append+" \""+PDBPATH+"\"; \
center; \
select */"+pdf+"; cpk off;wireframe off;delete solvent;"

    if(hasAX) {			// calculate AXis rotation only if present
        loadstr += "\
select */"+axf+"; cpk off; wireframe "+awf+"; color "+AXCOLOR+"; \
As={(*/"+axf+")[1]}.xyz; \
Af={(*/"+axf+")[0]}.xyz; \
df=As-Af; \
vf= @{df/sqrt(df*df)}; \
AyQ = !quaternion(cross({0 1 0},vf), acos({0 1 0}*vf));";
    }

    if(hasBB) {			// change representation of BB only if present
        loadstr += "\
select */"+bbf+"; cpk off; wireframe "+bwf+"; color gray;";
    }

    loadstr += "\
sps=off; \
piso=off; \
niso=off; \
frame all; \
print \"loaddone\"; \
";

    /*
     * create all GUI elements
     */
    jmolApplet([APPWIDTH, APPHEIGHT], loadstr);

    // show/hide models:
    //jmolHtml("<div id=\"isoload\"><img class=\"loader\" src=\""+LOADERIMG+"\" />Generating...</div>");
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
    return jmolScript("background "+el.id);
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
    
    if(msg.indexOf("loaddone") >= 0) { // done loading all models
	var chistr = jmolGetPropertyAsJSON("chaininfo"); // chainInfo object
	var chi = eval("("+chistr+")");
	_initPDB(chi.chaininfo);
	_nash(pdf, 1);
	_prsh(pdf, 1);
	if(hasBB) {    // only init BB if the backbone file is present
	    _initBB(chi.chaininfo);
	}
    } else if(msg.indexOf("isocreate") >= 0) { // start creating isosurface
	document.getElementById("isoload").style.display = "block";
    }else if(msg.indexOf("isodone") >= 0) { // done creating isosurface
	document.getElementById("isoload").style.display = "none";
    }
    return true;
}

//---------------
function _initPDB(chi) { // initialise dynamically-derived pdb: nucleic acid [and protein]

    var nacontrols = "";
    var procontrols = "";
    jmolSetDocument(null);
    // Nucleic acid
    nacontrols += jmolHtml("<div class=\"disptk\"><div class=\"nacheck\">");
    nacontrols += jmolCheckbox("javascript _nash("+pdf+",2)","javascript  _nash("+pdf+",0)", "Nucleic Acid: ", "checked");
    nacontrols += jmolHtml("</div><div class=\"naradio\">");
    pradio = [];
    for(var prep=0; prep<NDISPMODES.length; prep++) {
	var checked = 0+(prep==(0+(!hasBB)));
	pradio.push(["javascript _nash("+pdf+",1)", NDISPNAMES[prep], checked]);
    }
    nacontrols += jmolRadioGroup(pradio, null, "jmolRadioGroup"+naid+"disp")
    nacontrols += jmolHtml("</div></div>");

    // determine whether a protein is present
    if(jmolGetPropertyAsArray("atomList","protein AND */"+pdf).length > 0) {
	hasProtein=true;
	// Protein
	procontrols += jmolHtml("<div class=\"disptk\"><div class=\"nacheck\">");
	procontrols += jmolCheckbox("javascript _prsh("+pdf+",2)","javascript  _prsh("+pdf+",0)", "Protein: ", "checked");
	procontrols += jmolHtml("</div><div class=\"naradio\">");
	pradio = [];
	for(var prep=0; prep<PDISPMODES.length; prep++) {
	    var checked = 0;
	    pradio.push(["javascript _prsh("+pdf+",1)", PDISPNAMES[prep], checked]);
	}
	procontrols += jmolRadioGroup(pradio, null, "jmolRadioGroup"+prid+"disp")
	procontrols += jmolHtml("</div></div>");
    }else{
	hasProtein=false;
    }
	
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
    
    var bbstr = "";
    for(gi=0; gi<grooves.length; gi++) {
	bbstr += "select */"+bbf+" AND "+grooves[gi]+"; color bonds "+GRVCOLORS[gi]+";";
    }
    bbstr += "select */"+bbf+" AND "+bbsele+"; color "+BACCOLOR+";";
    var ret1 = jmolScript(bbstr);
    
    /*
     * Insert the controls for bb and grooves
     */
    
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
    document.getElementById("bbdisplay").innerHTML = bbcontrols;
    return ret1;
}
