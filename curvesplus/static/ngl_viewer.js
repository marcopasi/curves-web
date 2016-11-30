/*
 * Copyright (C) 2015-2016 Marco Pasi <mf.pasi@gmail.com> 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 */

/**** NGL for Curves+
  * version 0.3
  ****/

/**** Changelog
 * v0.1		First draft
 * v0.2		Reproduce 90% capabilities of the Jmol viewer
 * v0.3		RepresentationGroups API simplify GUI creation
 ****/

/**** TODO
 * - variable groove number
 * - align and spin wrt Axis
 * - reset orientation to initial
 */


/* GLOBAL SETTINGS */
var AWF = 0.3,                  // radius for Axis
    BWF = 0.3,                  // radius for Backbone
    GWF = 0.3,                  // radius for Grooves
    CWF = 0.3;                  // radius for Curvature vectors
var BGCOLORS =["lightgray","white","black","powderblue"];

var DEBUG = false;

/* GLOBALS */
var stage, repdata;

/* Create the viewer */
function ngl_viewer(AXPATH, BBPATH, CRPATH, PDBPATH) {
    repdata = {};
    stage = new NGL.Stage("viewport",
                          {"cameraType": "orthographic"});

    // Create RepresentationGroups for the input PDB
    var pdbRG = stage.loadFile(PDBPATH)
        .then(function(c) {return c.centerView();})
        .then(do_input, error);

    var axRG, bbRG, crRG;
    if(typeof(AXPATH) != "undefined") {
        // Create RepresentationGroups for the axis PDB
        axRG = stage.loadFile(AXPATH)
            .then(function(c) {return c.centerView();})
            .then(do_ax, error);
    }
    if(typeof(BBPATH) != "undefined") {
        // Create RepresentationGroups for the backbone PDB
        bbRG = stage.loadFile(BBPATH).then(do_bb, error);
    }
    if(typeof(CRPATH) != "undefined") {
        // Create RepresentationGroups for the curvature PDB
        crRG = stage.loadFile(CRPATH).then(do_cr, error);
    }

    // Wall: resolve all RepresentationGroups
    Promise.all([pdbRG, axRG, bbRG, crRG]).then(function(RG) {
        // Aggregate RepresentationGroups in repdata
        RG.forEach(function(rep) {$.extend(repdata, rep);});
        return repdata;
    }).then(
        // Write GUI for RepresentationGroups
        // in specific containers, in a specific order.
        function(RGdata) {
            console.log(RGdata);
            var lc = $("#"+"lcontrols");
            lc.append(RGdata["Nucleic Acid"].GUI("nadisplay"),
                      RGdata["Axis"].GUI("axdisplay"),
                      RGdata["Backbone"].GUI("bbdisplay"),
                      RGdata["Groove12"].GUI("gr1display"),
                      RGdata["Groove21"].GUI("gr2display"),
                      RGdata["Curvature"].GUI("crdisplay"));
            var rc = $("#"+"rcontrols");
            rc.append(RGdata["Protein"].GUI("prodisplay"));
            rc.append(GUI_extras());
        });
}

function GUI_extras() {
    /* Define extra GUI elements.
     */
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
    
    // Buttons
    var ddiv = $("<div/>", {"class": "buttons"});
    ddiv.append($("<input/>", {"type": "button",
                               "value": "Center"})
                .click(function(e) {stage.centerView();}),
                $("<br/>"),
                $("<input/>", {"type": "button",
                               "value": "Spin"})
                .click(function(e) {
                    if(stage.spin) {
                        stage.setSpin(null, 0);
                        stage.spin = undefined
                    } else {
                        stage.spin = true;
                        stage.setSpin([0,1,0], .01);
                    }}));

    return [cdiv, ddiv];
}

/* Representation callbacks
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
        new RepresentationGroup(comp, "Nucleic Acid", "nucleic",
                      [comp.addRepresentation( "licorice",   {"colorScheme": "uniform",
                                                              "colorValue":  "gray"}),
                       comp.addRepresentation( "ball+stick", {"colorScheme": "element"}),
                       comp.addRepresentation( "spacefill",  {"colorScheme": "uniform",
                                                              "colorValue":  "yellow"})]),
        // Protein
        "Protein":
        new RepresentationGroup(comp, "Protein", "protein",
                      [comp.addRepresentation( "cartoon", {"colorScheme":   "uniform",
                                                           "colorValue":    "steelblue"}),
                       comp.addRepresentation( "licorice", {"colorScheme":  "element"}),
                       comp.addRepresentation( "spacefill", {"colorScheme": "uniform",
                                                             "colorValue":  "deepskyblue"})])
    };
}

function do_ax(comp) {
    return {
        "Axis":
        new RepresentationGroup(comp, "Axis", null,
                                [comp.addRepresentation( "licorice", {"colorScheme": "uniform",
                                                                      "colorValue":  "blue",
                                                                      "radius":      AWF})])
    };
}

function do_bb(comp) {
    return {
        "Backbone": 
        new RepresentationGroup(comp, "Backbone", "(:A or :B)",
                      [comp.addRepresentation( "licorice", {"colorScheme": "uniform",
                                                            "colorValue":  "red",
                                                            "radius":      BWF})]),
        "Groove12": 
        new RepresentationGroup(comp, "Groove12", ":C",
                      [comp.addRepresentation( "licorice", {"colorScheme": "uniform",
                                                            "colorValue":  "lawngreen",
                                                            "radius":      GWF})]),
        "Groove21": 
        new RepresentationGroup(comp, "Groove21", ":D",
                      [comp.addRepresentation( "licorice", {"colorScheme": "uniform",
                                                            "colorValue":  "orange",
                                                            "radius":      GWF})]),
        
    };
}

function do_cr(comp) {
    return {
        "Curvature":
        new RepresentationGroup(comp, "Curvature", null,
                                [comp.addRepresentation( "licorice", {"colorScheme": "uniform",
                                                                      "colorValue":  "blue",
                                                                      "radius":      CWF})])
    };
}

/* Representation groups API
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
            self.addRepresentation(repr);
        });
};

RepresentationGroup.prototype.addRepresentation = function(repr) {
    // Apply default parameters
    repr.setParameters(this.defaultParameters);
    // Set selection if defined
    if(this.selection)
        repr.setSelection(this.selection);
    // Hide initially
    repr.setVisibility(false);
    this.reprList.push(repr);
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

RepresentationGroup.prototype.GUI = function(class_name) {
    /*
     * Write HTML to control the visibility of this group.
     *
     * "class_name" is the class used for the container DIV.
     */
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
                $("<label/>", {"for": radioname+"_"+i}).append(capitalize(cd.name)),
                "&nbsp;"
                ));
        });
        c.append(d);
    }
    this.update();
    return c;
}

/* Promise functions */
function error(err) {
    console.log(err);
}

/* Utilities */
function capitalize(string) {
    return string.charAt(0).toUpperCase() + string.slice(1);
}
