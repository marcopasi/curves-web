/*
 * Copyright (C) 2015-2016 Marco Pasi <mf.pasi@gmail.com> 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 */

/**** NGL for Curves+
  * version 0.9
  ****/

/**** Changelog
 * v0.1		First draft
 * v0.2		Reproduce 90% capabilities of the Jmol viewer
 * v0.3		RepresentationGroups API simplify GUI creation
 * v0.4		Reset orientation to initial
 * v0.5		Spin wrt/ Axis
 * v0.6		Align Axis vertical
 * v0.7		Screenshot button added
 * v0.8		color by nucleotide
 * v0.9		Debugged Safari download with hack (see TODO)
 * v0.9.1	Refactored color settings
 ****/

/**** TODO
 * - variable groove number
 * - remove Safari Hack for image download when Safari supports download file name setting (see webkit bug 102914).
 */


/* GLOBAL SETTINGS */
var AWF = 0.3,                  // radius for Axis
    BWF = 0.3,                  // radius for Backbone
    GWF = 0.3,                  // radius for Grooves
    CWF = 0.3;                  // radius for Curvature vectors
var Aco = "midnightblue",
    Bco = "darkred",
    Cco = Aco,
    Gcos= ["mediumvioletred","darkorange","pink","silver"],
    Ncw = "gray",
    Ncs = "khaki",
    Nso = 0.5, // opacity
    Pcc = "steelblue",
    Pcs = "skyblue",
    Pso = 0.5; // opacity
var BGCOLORS =["lightgray","white","black","powderblue"];
// NDB colors
var lowsaturation = //RGBYC
    [0xE6BEBE, 0xBEE6BE, 0xBEBEE6, 0xE6E6AA, 0xBEE6E6];
var higsaturation = //RGBYC
    [0xFF7070, 0xA0FFA0, 0xA0A0FF, 0xFFFF70, 0x70FFFF];
var rgbyc = higsaturation;
var NDBColors = NGL.ColorMakerRegistry.addSelectionScheme( [ // A red, T blue, C yellow, G green, and U cyan. 
    [rgbyc[0],"DA or A"],
    [rgbyc[1],"DG or G"],
    [rgbyc[2],"DT"],
    [rgbyc[3],"DC or C"],
    [rgbyc[4],"U"],
    ["gray","*"]
]);

var DEBUG = false;

/* GLOBALS */
var stage, repdata, dna_axis, orientation, zoom;


/*************************
 * Create the viewer
 */
function ngl_viewer(AXPATH, BBPATH, CRPATH, PDBPATH) {
    repdata = {};
    stage = new NGL.Stage("viewport",
                          {"cameraType": "orthographic",
                           "backgroundColor": "white"});

    // Create RepresentationGroups for the input PDB
    var pdbRG = stage.loadFile(PDBPATH)
        .then(do_input, error);

    var axRG, bbRG, crRG;
    // Define dummy axis if we lack one
    dna_axis = new NGL.Vector3(0,1,0);
    if(AXPATH != "") {
        // Create RepresentationGroups for the axis PDB
        axRG = stage.loadFile(AXPATH)
            .then(function(c) {
                // Get Axis approximate axis
                dna_axis = get_axis(c.structure);
                return c.centerView();})
            .then(do_ax, error);
    } else {
        pdbRG.then(function(rg) {rg["Nucleic Acid"].component.centerView(); return rg;});
    }
    if(BBPATH != "") {
        // Create RepresentationGroups for the backbone PDB
        bbRG = stage.loadFile(BBPATH).then(do_bb, error);
    }
    if(CRPATH != "") {
        // Create RepresentationGroups for the curvature PDB
        crRG = stage.loadFile(CRPATH).then(do_cr, error);
    }

    // Wall: resolve all RepresentationGroups
    Promise.all([pdbRG, axRG, bbRG, crRG]).then(function(RG) {
        // Get initial orientation and zoom
        orientation = stage.viewer.getOrientation();
        zoom = stage.viewer.camera.zoom;
        // Aggregate RepresentationGroups in repdata
        RG.forEach(function(rep) {$.extend(repdata, rep);});
        return repdata;
    }).then(
        // Write GUI for RepresentationGroups
        // in specific containers, in a specific order.
        function(RGdata) {
            var lc = $("#"+"lcontrols");
            if(RGdata["Nucleic Acid"])
                lc.append(RGdata["Nucleic Acid"].GUI("nadisplay"));
            if(RGdata["Axis"])
                lc.append(RGdata["Axis"].GUI("axdisplay"));
            if(RGdata["Backbone"])
                lc.append(RGdata["Backbone"].GUI("bbdisplay"));
            if(RGdata["Groove12"])
                lc.append(RGdata["Groove12"].GUI("gr1display"));
            if(RGdata["Groove21"])
                lc.append(RGdata["Groove21"].GUI("gr2display"));
            if(RGdata["Curvature"])
                lc.append(RGdata["Curvature"].GUI("crdisplay"));
            
            var rc = $("#"+"rcontrols");
            if(RGdata["Protein"])
                rc.append(RGdata["Protein"].GUI("prodisplay"));
            rc.append(GUI_extras(RGdata["Axis"]));
            more_GUI_extras();
        });

    window.addEventListener(
        "resize", function( event ){
            stage.handleResize();
        }, false
    );
}

/*************************
 * Define extra GUI elements.
 */
function GUI_extras(axis) {
    // Background
    var cdiv = $("<div/>", {"class": "colors"});
    cdiv.append("Background: ");
    BGCOLORS.forEach(function(c) {
        cdiv.append(
            $("<a/>", {"class": "dummylink"}).append(
                $("<div/>", {"id": "",
                             "class": "cimg"})
                    .css("background-color", c)
                    .click(function(e) {stage.viewer.setBackground(c);})));
    });

    if(!axis) return cdiv;
    
    function spin() {
        // spin stage around DNA axis
        if(stage.spin) {
            stage.viewer.setSpin(undefined, 0);
            stage.spin = undefined
        } else {
            stage.viewer.setSpin(undefined, Math.PI/120);
            stage.spin = true;
        }
    }

    stage.viewer.setSpin(undefined, 0); // set no spinning
    // Keep spinning axis aligned with DNA axis
    function setspin() {
        stage.viewer.setSpin(invcc(dna_axis), undefined);
    }
    stage.viewer.controls.addEventListener("change", setspin);
    setspin();                  // set axis first time
    
    function align() {
        //XXX aligns axis horizontal
        rotateCameraTo(dna_axis);
        rotateCameraAxisAngle(cc(new NGL.Vector3(1,0,0)), -Math.PI/2)
    }

    // requires previously defined initial orientation and zoom
    function orient() {
        // reset orientation to initial
        stage.viewer.camera.zoom = zoom;
        stage.viewer.setOrientation(orientation);
    }
    
    // Buttons
    var ddiv = $("<div/>", {"class": "buttons"});
    ddiv.append(
        $("<input/>", {"type": "button",
                       "value": "Reset orientation"})
            .click(orient),
        $("<br/>"),
        $("<input/>", {"type": "button",
                       "value": "Align Axis"})
            .click(align),
        $("<br/>"),
        $("<input/>", {"type": "button",
                       "value": "Spin"})
            .click(spin));

    return [cdiv, ddiv];
}

function safariw(data, target) {
    var url = URL.createObjectURL( data );
    target.location.href = url;
}

function more_GUI_extras() {
    // Safari Hack: open image in new window
    if( typeof window === "undefined" ) return false;
    var ua = window.navigator.userAgent,
        isSafari = ( /Safari/i.test( ua ) && ! /Chrome/i.test( ua ) ),
        safariwin = null;
    console.log(isSafari);
    // More Buttons
    function image() {
        if(isSafari) { // Safari Hack: open window early, otherwise Safari will block it!
            safariwin = window.open();
            safariwin.document.write("<html><head><link rel='stylesheet' href='/static/style.css'></head><body><h1>Curves+ web server</h1><h3>Thank you for your patience...</h3>Please wait while screenshot is being created. This may take a few seconds...</body></html>")
        }
        var fname = "screenshot.png"
        stage.makeImage({
            factor: 4,
            antialias: true,
            trim: false,
            transparent: false
        }).then( function( blob ){
            if(isSafari) { // Safari Hack: set new window URL to image
                return safariw(blob, safariwin);
            }
            return NGL.download( blob, fname );
        });
    }

    var ediv = $("<div/>", {"style": "position:absolute; bottom: 5px; right: 5px;"});
    ediv.append(
        $("<img/>", {"src": "/static/img/camera.svg",
                     "width": "25px"})
            .click(image));
    $(stage.viewer.container)
        .css("position", "relative")
        .append(ediv);
}

/*************************
 * Representation callbacks
 *
 * Configure representations here.
 * Each method creates a dictionary of RepresentationGroups, one
 * for each selection relevant for the specified component.
 * These are subsequently aggregated and used to design GUI.
 *
 */   
function do_input(comp) {
    return {
        // Nucleic
        "Nucleic Acid":
        new RepresentationGroup(comp, "Nucleic Acid", "nucleic")
            .addRepresentation( "Gray",
                                comp.addRepresentation( "licorice",   {"colorScheme": "uniform",
                                                                       "colorValue":  Ncw}))
            .addRepresentation( "Color",
                                comp.addRepresentation( "licorice",   {"colorScheme": NDBColors}))
            .addRepresentation( "Ball&Stick",
                                comp.addRepresentation( "ball+stick", {"colorScheme": "element"}))
            .addRepresentation( "Surface",
                                comp.addRepresentation( "surface",    {"opacity": Nso,
                                                                       "colorScheme": "uniform",
                                                                       "colorValue":  Ncs})),
        // Protein
        "Protein":
        new RepresentationGroup(comp, "Protein", "protein")
            .addRepresentation( "Cartoon",
                                comp.addRepresentation( "cartoon",  {"colorScheme":   "uniform",
                                                                    "colorValue":    Pcc}))
            .addRepresentation( "Wire",
                                comp.addRepresentation( "licorice", {"colorScheme":  "element"}))
            .addRepresentation( "Surface",
                                comp.addRepresentation( "surface",  {"opacity": Pso,
                                                                     "colorScheme": "uniform",
                                                                     "colorValue":  Pcs}))
    };
}

function do_ax(comp) {
    return {
        "Axis":
        new RepresentationGroup(comp, "Axis", null)
            .addRepresentation( null,
                               comp.addRepresentation( "licorice", {"colorScheme": "uniform",
                                                                    "colorValue":  Aco,
                                                                    "radius":      AWF}))
    };
}

function do_bb(comp) {
    return {
        "Backbone": 
        new RepresentationGroup(comp, "Backbone", "(:A or :B)")
            .addRepresentation(null,
                               comp.addRepresentation( "licorice", {"colorScheme": "uniform",
                                                                    "colorValue":  Bco,
                                                                    "radius":      BWF})),
        "Groove12": 
        new RepresentationGroup(comp, "Groove12", ":C")
            .addRepresentation(null,
                               comp.addRepresentation( "licorice", {"colorScheme": "uniform",
                                                                    "colorValue":  Gcos[0],
                                                                    "radius":      GWF})),
        "Groove21": 
        new RepresentationGroup(comp, "Groove21", ":D")
            .addRepresentation(null,
                               comp.addRepresentation( "licorice", {"colorScheme": "uniform",
                                                                    "colorValue":  Gcos[1],
                                                                    "radius":      GWF})),
    };
}

function do_cr(comp) {
    return {
        "Curvature":
        new RepresentationGroup(comp, "Curvature", null)
            .addRepresentation(null,
                               comp.addRepresentation( "licorice", {"colorScheme": "uniform",
                                                                    "colorValue":  Cco,
                                                                    "radius":      CWF}))
    };
}

/*************************
 * Representation groups API
 */
RepresentationGroup = function(component, name, selection = null, representations = null, 
                               defaultParameters = {}) {
    
    /* Representation Group
     * 
     * Helps define groups of mutually-exclusive representations, and
     * enables writing simple HTML GUI to control visibility.
     * Uses the "enabled" property of NGL.Representation to keep track
     * of user-specified visibility.
     *
     * Arguments
     * ---------
     *
     * component	 The NGL.Component of the representations
     * name		 A string describing the group (used in GUI)
     * representations	 Optional list of representations to add
     * selection	 A selection of the subset of component atoms
     * 			 to which the representation is applied.
     * defaultParameters A dictionary of NGL.Representation parameters
     *			 to be applied to all representations.
     */
    this.component = component;
    this.name = name;
    this.selection = selection;
    this.defaultParameters = defaultParameters;
    this.enabled = true;
    
    this.reprList = [];

    var self = this;
    if(representations) 
        representations.forEach(function(repr) {
            self.addRepresentation(null, repr);
        });
};

RepresentationGroup.prototype.addRepresentation = function(display_name, repr) {
    repr.display_name = display_name;
    // Apply default parameters
    repr.setParameters(this.defaultParameters);
    // Set selection if defined
    if(this.selection)
        repr.setSelection(this.selection);
    // Hide initially
    repr.setVisibility(false);
    this.reprList.push(repr);
    return this;
}
    
RepresentationGroup.prototype.toggle = function(checked) {
    // Toggle the enabled state of this group
    this.enabled = checked;
    this.update();
}

RepresentationGroup.prototype.hideall = function() {
    // Hide all representations in group
    this.reprList.forEach(function(cd) {
        cd.setVisibility(false);
    });    
}

RepresentationGroup.prototype.enable = function(ci=-1) {
    // Enable one representation of the group, by index
    this.reprList.forEach(function(cd, i) {
        cd.enabled = i == ci;
    });
    this.update();
}

RepresentationGroup.prototype.update = function() {
    // Update representation visilibity
    if(!this.enabled) this.hideall();
    else {
        this.reprList.forEach(function(cd) {
            cd.setVisibility(cd.enabled);
        });
    }
}

RepresentationGroup.prototype.all_empty = function() {
    // Check if all representations in group are empty
    return this.reprList.every(function(repr) {
        return repr.repr.structureView.atomCount == 0;
    });
}

RepresentationGroup.prototype.GUI = function(class_name) {
    /*
     * Write HTML to control the visibility of this group.
     *
     * "class_name" is the class used for the container DIV.
     */
    
    // If all the representation are empty, return empty GUI
    if(this.all_empty()) return null;
    
    var self = this,
        c = $("<div/>", {"class": class_name})
    
    c.append($("<div/>").append(
        $("<input/>", {"type": "checkbox",
                       "id": class_name,
                       "checked": true})
            .click(function(e) {self.toggle(this.checked);}),
        $("<label/>", {"for": class_name}).append(this.name)
    ));
    
    // Enable first (only) representation
    this.reprList[0].enabled = true;
    
    if(this.reprList.length > 1) {
        var d = $("<div/>", {"class":"naradio"}),
            radioname = class_name+"radio";
        this.reprList.forEach(function(cd, i) {
            cd.enabled = (i == 0) ? true : false;
            d.append($("<span/>").append(
                $("<input/>", {"type": "radio",
                               "name": radioname,
                               "id": radioname+"_"+i,
                               "checked": i==0})
                    .click(function(e) {self.enable(i);}),
                $("<label/>", {"for": radioname+"_"+i}).append(capitalize(cd.display_name)),
                "&nbsp;"
                ));
        });
        c.append(d);
    }
    this.update();
    return c;
}

/*************************
 * Promise functions
 */
function error(err) {
    console.log(err);
}

/*************************
 * Utilities
 */
function capitalize(string) {
    return string.charAt(0).toUpperCase() + string.slice(1);
}

function get_axis(structure) {
    var atoms = structure.atomStore,
        n = atoms.count-1,
        x = atoms.x[0] - atoms.x[n],
        y = atoms.y[0] - atoms.y[n],
        z = atoms.z[0] - atoms.z[n];
    return new NGL.Vector3(x,y,z).normalize();
}

function rotateCameraTo(end) {
    var camera = stage.viewer.camera,
        start = camera.getWorldDirection(),
        target = stage.viewer.controls.target;

    var angle = Math.acos(start.dot(end) / start.length() / end.length()),
        raxis = (new NGL.Vector3()).crossVectors(start, end).normalize();
    
    return rotateCameraAxisAngle(raxis, angle, target);
}

function rotateCameraAxisAngle(axis, angle, target) {
    if (!angle) return;
    if (!target) target = stage.viewer.controls.target;
    
    var camera = stage.viewer.camera,
        _eye = new NGL.Vector3();

    var quaternion = (new NGL.Quaternion()).setFromAxisAngle(axis, angle);
    // rotate the distance vector (_eye)
    _eye.subVectors(camera.position, target);
    _eye.applyQuaternion(quaternion);
    // rotate the camera's up vector
    camera.up.applyQuaternion(quaternion);
    // re-apply rotated distance to camera
    camera.position.addVectors(_eye, target);
    // reorient camera towards target
    camera.lookAt(target);
}

function cc(axis) {
    var camera = stage.viewer.camera,
        controls = stage.viewer.controls;

    var eye = (new NGL.Vector3()).copy( camera.position ).sub( controls.target ),
        eyeDirection = (new NGL.Vector3()).copy( eye ).normalize(),
        upDirection = (new NGL.Vector3()).copy( camera.up ).normalize(),
        sidewaysDirection = (new NGL.Vector3()).crossVectors( upDirection, eyeDirection ).normalize(),
        moveDirection = new NGL.Vector3();
    
    // console.log("  s=[", sidewaysDirection.toArray().toString(),
    //             "]; u=[", upDirection.toArray().toString(),
    //             "]; e=[", eyeDirection.toArray().toString(),
    //             "];");

    /*
     * The following operations are equivalent to:
     * 
     * moveDirection = M*axis
     *
     * where M is defined as follows:
     */
    
    var M = (new NGL.Matrix4()).makeBasis(sidewaysDirection, upDirection.clone().negate(), eyeDirection);
    
    eyeDirection.setLength( axis.z );
    upDirection.setLength( axis.y );
    sidewaysDirection.setLength( axis.x );

    // console.log("S", sidewaysDirection,
    //             "U", upDirection,
    //             "E", eyeDirection);
    
    moveDirection.copy( sidewaysDirection.sub( upDirection ).add( eyeDirection ) );

    // console.log("D",moveDirection2axis.clone().applyMatrix4(M).normalize().sub(moveDirection));
    
    return moveDirection;
}

function invcc(axis) {
    // Implement inverse operation performed in NGL.Viewer.rotate
    // to directly define axis.
    var camera = stage.viewer.camera,
        controls = stage.viewer.controls;

    var eye = (new NGL.Vector3()).copy( camera.position ).sub( controls.target ),
        eyeDirection = (new NGL.Vector3()).copy( eye ).normalize(),
        upDirection = (new NGL.Vector3()).copy( camera.up ).normalize(),
        sidewaysDirection = (new NGL.Vector3()).crossVectors( upDirection, eyeDirection ).normalize();
    
    var M = (new NGL.Matrix4()).makeBasis(sidewaysDirection, upDirection.clone().negate(), eyeDirection);
    return axis.clone().applyMatrix4(M.transpose()).normalize();
}
