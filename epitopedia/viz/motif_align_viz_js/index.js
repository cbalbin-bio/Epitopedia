// feed from json. make json request from data in meta tag in html

var query = data["EPI_SEQ Input Structure"] + ".cif";
var query_nums = data["SeqBMM Input Struc Res Nums"];
var mimic_nums = data["EPI_PDB Rep Res Nums"].map(Number);
var mimic = data["EPI_PDB Rep PDB"] + ".cif";
var aln = data["EPI_PDB TMalign PDB"];

var query_chain = query.replace(".cif", "").split("_").at(-1);
var mimic_chain = mimic.replace(".cif", "").split("_").at(-1);

var query_selection = `${query_nums.at(0)}-${query_nums.at(-1)}:${query_chain}`
var mimic_selection = `${mimic_nums.at(0)}-${mimic_nums.at(-1)}:${mimic_chain}`

var mimic_plddt = data["mmCIF_SEQ lplddt"].split(" ").map(Number)
var afScale = chroma.scale(['orange', 'yellow', 'cyan', 'blue']).domain([0.0, 100.0]);


function colorSchemeAln(af, query_selected, mimic_selected) {
    var nglColorScheme = NGL.ColormakerRegistry.addScheme(function (params) {
        this.atomColor = function (atom) {
            rindex = atom.resno;

            // console.log(atom.chainname)



            if (atom.chainname == "B") {
                if (af) {
                    if (rindex == mimic_selected) {
                        return "0x0000ff"
                    } else {
                        return afScale(mimic_plddt[atom.residueIndex]).hex().replace("#", "0x");
                    }
                }
                else {
                    if (rindex == mimic_selected) {
                        return "0x0000ff"
                    } else {
                        return "0x7899d2"
                    }

                }
            }
            else if (atom.chainname == "A") {
                if (rindex == query_selected) {
                    return "0x0000ff"
                } else {
                    return "0xcc955d"
                }

            }
        }
    })




    return nglColorScheme
};



function colorScheme(start, end, chain, type, af, selected) {
    var nglColorScheme = NGL.ColormakerRegistry.addScheme(function (params) {
        this.atomColor = function (atom) {
            rindex = atom.resno;

            // console.log(atom.chainname)
            if (af) {
                if (rindex == selected && atom.chainname == chain) {
                    return "0x0000ff"
                } else if (rindex >= start && rindex <= end && atom.chainname == chain) {
                    return "0x00ff00";
                } else {
                    return afScale(mimic_plddt[atom.residueIndex]).hex().replace("#", "0x");
                }

            } else {
                if (type == "mimic") {
                    if (rindex == selected && atom.chainname == chain) {
                        return "0x0000ff"
                    } else if (rindex >= start && rindex <= end && atom.chainname == chain) {
                        return "0x00ff00";
                    } else {
                        return "0x7899d2"
                    }


                } else if (type == "input") {
                    if (rindex == selected && atom.chainname == chain) {
                        return "0x0000ff"
                    } else if (rindex >= start && rindex <= end && atom.chainname == chain) {
                        return "0x00ff00";
                    } else {
                        return "0xcc955d"

                    }
                }





            };
        };
    });
    return nglColorScheme
};






var schemeIdplddt = NGL.ColormakerRegistry.addScheme(function (params) {
    this.atomColor = function (atom) {
        rindex = atom.residueIndex;

        return scale(mimic_plddt[rindex]).hex().replace("#", "0x");

    };
});



var schemeIdQuery = NGL.ColormakerRegistry.addSelectionScheme([
    ["red", query_selection],
    ["#cc955d", `:${query_chain}`]
], "motif");


var schemeIdMimic = NGL.ColormakerRegistry.addSelectionScheme([
    ["red", mimic_selection],
    ["#7899d2", `:${mimic_chain}`]
], "motif");

var schemeIdAln = NGL.ColormakerRegistry.addSelectionScheme([
    ["#cc955d", ":A"],
    ["#7899d2", ":B"]
], "motif");







var stage1 = new NGL.Stage("viewport1", { backgroundColor: "black" });
stage1.loadFile("cif/" + query).then(function (o) {
    o.addRepresentation("cartoon", { color: colorScheme(query_nums.at(0), query_nums.at(-1), query_chain, "input", false, false) });
    o.addRepresentation("ball+stick", { sele: query_selection, color: colorScheme(query_nums.at(0), query_nums.at(-1), query_chain, "input", false, false) });
    o.autoView(query_selection);
    //o.autoView("22-26:A");

})
// stage1.loadFile("6XR8.cif",).then(function (o) {
//     o.addRepresentation("cartoon", { color: schemeId });
//     o.addRepresentation("surface", { sele: "22-26", opacity: 0.6, color: schemeId }); // pass schemeId here

//     var ap = o.structure.getAtomProxy()
//     var shape = new NGL.Shape("shape");
//     ap.residueIndex = 22 + Math.floor((26 - 22) / 2);



//     shape.addArrow([ap.positionToVector3()["x"], ap.positionToVector3()["y"], ap.positionToVector3()["z"]], [ap.positionToVector3()["x"], ap.positionToVector3()["y"], ap.positionToVector3()["z"] + 5], [1, 0, 1], 1.0)
//     var shapeComp = stage1.addComponentFromObject(shape);
//     console.log(stage1)
//     shapeComp.addRepresentation("buffer");

//     o.autoView();

// });

// stage2.compList[0].reprList[2].setSelection(2-20).setColor("red")
var stage2 = new NGL.Stage("viewport2", { backgroundColor: "black" });
stage2.loadFile("cif/" + mimic).then(function (o) {

    o.addRepresentation("cartoon", { color: colorScheme(mimic_nums.at(0), mimic_nums.at(-1), mimic_chain, "mimic", data["mmCIF_SEQ AF"], false), })
    o.addRepresentation("ball+stick", { sele: mimic_selection, color: colorScheme(mimic_nums.at(0), mimic_nums.at(-1), mimic_chain, "mimic", data["mmCIF_SEQ AF"], false) });

    // o.addRepresentation("cartoon", { sele: mimic_selection, color: colorScheme(mimic_nums.at(0), mimic_nums.at(-1), mimic_chain, "input", true, false) });

    // o.reprList.at(-2).toggleVisibility()


    o.autoView(mimic_selection);
});


// var stage3 = new NGL.Stage("viewport3", { backgroundColor: "white" });
// stage3.loadFile("6XR8.cif",).then(function (o) {
//     o.addRepresentation("validation", { color: schemeIdQuery });  // pass schemeId here
//     o.autoView();
// });


var stage3 = new NGL.Stage("viewport3", { backgroundColor: "black" });
stage3.loadFile("aln/" + aln).then(function (o) {
    o.addRepresentation("ball+stick", { color: colorSchemeAln(data["mmCIF_SEQ AF"], false) });  // pass schemeId here
    o.autoView();
});


window.addEventListener("resize", function (event) {
    stage1.handleResize();
    stage2.handleResize();
    // stage3.handleResize();
    stage3.handleResize();
}, false);


var elms = document.getElementsByClassName('residue');
for (const elm of elms) {
    elm.addEventListener("mouseenter", function (event) {

        // highlight the mouseleave target
        event.target.getElementsByTagName("rect")[0].style.fill = "blue";

        for (const repr of stage1.compList[0].reprList) {
            repr.setColor(
                colorScheme(query_nums.at(0), query_nums.at(-1), event.target.dataset.chainQuery, "input", false, event.target.dataset.resnumQuery));
            repr.update({ color: true });
        };

        for (const repr of stage2.compList[0].reprList) {
            repr.setColor(
                colorScheme(mimic_nums.at(0), mimic_nums.at(-1), event.target.dataset.chainMimic, "mimic", data["mmCIF_SEQ AF"], event.target.dataset.resnumMimic));
            repr.update({ color: true });
        };


        for (const comp of stage3.compList) {
            comp.reprList[0].setColor(colorSchemeAln(data["mmCIF_SEQ AF"], event.target.dataset.resnumQuery, event.target.dataset.resnumMimic)


            );
            comp.reprList[0].update({ color: true });
        };


        // reset the color after a short delay



        // console.log(elm.dataset);

    }, false);


    elm.addEventListener("mouseleave", function (event) {
        // highlight the mouseleave target
        event.target.getElementsByTagName("rect")[0].style.fill = "green";

        // reset the color after a short delay
        for (const repr of stage1.compList[0].reprList) {
            repr.setColor(
                colorScheme(query_nums.at(0), query_nums.at(-1), event.target.dataset.chainQuery, "input", false, false)
            );
            repr.update({ color: true });
        };


        for (const repr of stage2.compList[0].reprList) {
            repr.setColor(
                colorScheme(mimic_nums.at(0), mimic_nums.at(-1), event.target.dataset.chainMimic, "mimic", data["mmCIF_SEQ AF"], false))

            repr.update({ color: true });
        };


        for (const repr of stage3.compList[0].reprList) {
            repr.setColor(colorSchemeAln(data["mmCIF_SEQ AF"], false))



            repr.update({ color: true });
        };



        // console.log(elm.dataset);

    }, false);

};
