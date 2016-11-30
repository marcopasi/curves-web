/*
 * Copyright (C) 2015-2016 Marco Pasi <mf.pasi@gmail.com> 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 */

/**** NGL for Curves+
  * version 0.2
  ****/

/**** TODO
 * - variable groove number
 * - 
 * - 
 * - 
 */

/* GLOBAL SETTINGS */
var AWF = 0.3,                  // radius for Axis
    BWF = 0.3,                  // radius for Backbone
    GWF = 0.3,                  // radius for Grooves
    CWF = 0.3;                  // radius for Curvature vectors
var BGCOLORS =["lightgray","white","black","powderblue"];

var DEBUG = false;

/* GLOBALS */
var stage, compdata;

/* Create the viewer */
function ngl_viewer(AXPATH, BBPATH, CRPATH, PDBPATH) {
    compdata = {};
    stage = new NGL.Stage("viewport",
                          {"cameraType": "orthographic"});
    
    var pdbcomp = stage.loadFile(PDBPATH)
        .then(function(c) {return c.centerView();})
        .then(do_input, error);

    var axcomp, bbcomp, crcomp;
    if(typeof(AXPATH) != "undefined") {
        axcomp = stage.loadFile(AXPATH)
            .then(function(c) {return c.centerView();})
            .then(do_ax, error);
    }
    if(typeof(BBPATH) != "undefined") {
        bbcomp = stage.loadFile(BBPATH).then(do_bb, error);
    }
    if(typeof(CRPATH) != "undefined") {
        crcomp = stage.loadFile(CRPATH).then(do_cr, error);
    }
    
    Promise.all([pdbcomp, axcomp, bbcomp, crcomp]).then(function(cdata) {
        cdata.forEach(function(cd) {$.extend(compdata, cd);});
        return compdata
    }).then(
        function(cdata) {
            var lc = $("#"+"lcontrols");
            lc.append(r_GUI(cdata["Nucleic Acid"], "Nucleic Acid", "nadisplay"),
                      r_GUI(cdata["Axis"], "Axis", "axdisplay"),
                      r_GUI(cdata["Backbone"], "Backbone", "bbdisplay"),
                      r_GUI(cdata["Groove12"], "Groove12", "gr1display"),
                      r_GUI(cdata["Groove21"], "Groove21", "gr2display"),
                      r_GUI(cdata["Curvature"], "Curvature", "crdisplay"));
            var rc = $("#"+"rcontrols");
            rc.append(r_GUI(cdata["Protein"], "Protein", "prodisplay"));
            rc.append(GUI_extras());
        });
}

function GUI_extras() {
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
 */
function do_input(comp) {
    return {
        // Nucleic
        "Nucleic Acid":
        r_set_selection("nucleic",
                      [comp.addRepresentation( "licorice",   {"colorScheme":"uniform",
                                                              "colorValue":"gray"}),
                       comp.addRepresentation( "ball+stick", {"colorScheme":"element"}),
                       comp.addRepresentation( "spacefill",  {"colorScheme":"uniform",
                                                              "colorValue":"yellow"})]),
        // Protein
        "Protein":
        r_set_selection("protein",
                      [comp.addRepresentation( "cartoon", {"colorScheme":"uniform",
                                                           "colorValue":"steelblue"}),
                       comp.addRepresentation( "licorice", {"colorScheme":"element"}),
                       comp.addRepresentation( "spacefill", {"colorScheme":"uniform",
                                                             "colorValue":"deepskyblue"})])
    };
}

function do_ax(comp) {
    return {
        "Axis":
        [comp.addRepresentation( "licorice", {"colorScheme":"uniform",
                                              "colorValue":"blue",
                                              "radius": AWF})]
    };
}

function do_bb(comp) {
    return {
        "Backbone": 
        r_set_selection("(:A or :B)",
                      [comp.addRepresentation( "licorice", {"colorScheme":"uniform",
                                                            "colorValue":"red",
                                                            "radius": BWF})]),
        "Groove12": 
        r_set_selection(":C",
                      [comp.addRepresentation( "licorice", {"colorScheme":"uniform",
                                                            "colorValue":"lawngreen",
                                                            "radius": GWF})]),
        "Groove21": 
        r_set_selection(":D",
                      [comp.addRepresentation( "licorice", {"colorScheme":"uniform",
                                                            "colorValue":"orange",
                                                            "radius": GWF})]),
        
    };
}

function do_cr(comp) {
    return {
        "Curvature":
        [comp.addRepresentation( "licorice", {"colorScheme":"uniform",
                                              "colorValue":"blue",
                                              "radius": CWF})]
    };
}

/* Representation groups API
 *
 * TODO refactor in a class
 *  - the constructor could receive the component and optional selection
 *  - adding Representations should trigger:
 *    - apply default parameters
 *    - set_selection
 *    - making invisible
 */
function r_toggle(cdata, checked) {
    if(checked) r_update(cdata);
    else r_hideall(cdata);
}

function r_hideall(cdata) {
    cdata.forEach(function(cd) {
        cd.setVisibility(false);
    });    
}

function r_enable(cdata, ci=-1) {
    cdata.forEach(function(cd, i) {
        cd.enabled = i == ci;
    });
    r_update(cdata);
}

function r_update(cdata) {
    cdata.forEach(function(cd) {
        cd.setVisibility(cd.enabled);
    });
}

function r_set_selection(sele, comps) {
    comps.forEach(function(c) {c.setSelection(sele);});
    return comps;
}

function r_GUI(cdata, label, name) {
    /*
     *
     */
    if(cdata == undefined) return;
    var c = $("<div/>", {"class": name});
    c.append($("<div/>").append(
        $("<input/>", {"type": "checkbox",
                       "id": name,
                       "checked": true})
            .click(function(e) {r_toggle(cdata, this.checked);}),
        $("<label/>", {"for": name}).append(label)
    ));
    
    // Enable first (only) representation
    cdata[0].enabled = true;
    
    if(cdata.length > 1) {
        var d = $("<div/>", {"class":"naradio"}),
            radioname = name+"radio";
        cdata.forEach(function(cd, i) {
            cd.enabled = (i == 0) ? true : false;
            d.append($("<span/>").append(
                $("<input/>", {"type": "radio",
                               "name": radioname,
                               "id": radioname+"_"+i,
                               "checked": i==0})
                    .click(function(e) {r_enable(cdata, i);}),
                $("<label/>", {"for": radioname+"_"+i}).append(capitalize(cd.name)),
                "&nbsp;"
                ));
        });
        c.append(d);
    }
    r_update(cdata);
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
