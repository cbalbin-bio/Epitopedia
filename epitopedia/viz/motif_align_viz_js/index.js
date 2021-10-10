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
var scale = chroma.scale(['orange', 'yellow', 'cyan', 'blue']).domain([0.0, 100.0]);


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
    o.addRepresentation("cartoon", { color: schemeIdQuery });
    o.addRepresentation("ball+stick", { sele: query_selection, color: schemeIdQuery });
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
    if (data["mmCIF_SEQ AF"]) {
        o.addRepresentation("cartoon", { color: schemeIdplddt, })
    } else {
        o.addRepresentation("cartoon", { color: schemeIdMimic, });
    }
    o.addRepresentation("cartoon", { sele: mimic_selection, color: "red" });
    o.addRepresentation("ball+stick", { sele: mimic_selection, color: "red" });
    o.reprList.at(-2).toggleVisibility()


    o.autoView(mimic_selection);
});


// var stage3 = new NGL.Stage("viewport3", { backgroundColor: "white" });
// stage3.loadFile("6XR8.cif",).then(function (o) {
//     o.addRepresentation("validation", { color: schemeIdQuery });  // pass schemeId here
//     o.autoView();
// });


var stage3 = new NGL.Stage("viewport3", { backgroundColor: "black" });
stage3.loadFile("aln/" + aln).then(function (o) {
    o.addRepresentation("ball+stick", { color: schemeIdAln });  // pass schemeId here
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
        event.target.getElementsByTagName("rect")[0].style.fill = "green";

        for (const repr of stage1.compList[0].reprList) {
            repr.setColor(NGL.ColormakerRegistry.addSelectionScheme([
                ["green", event.target.dataset.resnumQuery + ":" + event.target.dataset.chainQuery],
                ["red", query_selection],

                ["#cc955d", `:${query_chain}`]
            ], "highlight"));
            repr.update({ color: true });
        };

        for (const repr of stage2.compList[0].reprList) {
            repr.setColor(NGL.ColormakerRegistry.addSelectionScheme([
                ["green", event.target.dataset.resnumMimic + ":" + event.target.dataset.chainMimic],
                ["red", mimic_selection],

                [schemeIdplddt, `:${mimic_chain}`]
            ], "highlight"));
            repr.update({ color: true });
        };


        for (const comp of stage3.compList) {
            comp.reprList[0].setColor(NGL.ColormakerRegistry.addSelectionScheme([
                ["green", event.target.dataset.resnumQuery + ":A"],
                ["green", event.target.dataset.resnumMimic + ":B"],
                ["#cc955d", ":A"],
                ["#7899d2", ":B"]
            ], "highlight"));
            comp.reprList[0].update({ color: true });
        };


        // reset the color after a short delay



        // console.log(elm.dataset);

    }, false);


    elm.addEventListener("mouseleave", function (event) {
        // highlight the mouseleave target
        event.target.getElementsByTagName("rect")[0].style.fill = "blue";

        // reset the color after a short delay
        for (const repr of stage1.compList[0].reprList) {
            repr.setColor(NGL.ColormakerRegistry.addSelectionScheme([
                ["red", event.target.dataset.resnumQuery + ":" + event.target.dataset.chainQuery],
                ["red", query_selection],

                ["#cc955d", `:${query_chain}`]

            ], "highlight"));
            repr.update({ color: true });
        };


        for (const repr of stage2.compList[0].reprList) {
            repr.setColor(NGL.ColormakerRegistry.addSelectionScheme([
                ["red", event.target.dataset.resnumMimic + ":" + event.target.dataset.chainMimic],
                ["red", mimic_selection],

                [schemeIdplddt, `:${mimic_chain}`]
            ], "highlight"));
            repr.update({ color: true });
        };


        for (const repr of stage3.compList[0].reprList) {
            repr.setColor(NGL.ColormakerRegistry.addSelectionScheme([
                ["#cc955d", ":A"],
                ["#7899d2", ":B"]
            ], "highlight"));
            repr.update({ color: true });
        };



        // console.log(elm.dataset);

    }, false);

};
