from flask import Flask, render_template, make_response, send_file
from jinja2 import Environment, FileSystemLoader
from rich import print

from epitopedia.app import config
from gemmi import cif
import sys
import json
import os


# def write_html(output_path, data):
#     with open(output_path, "w") as handle:

#         env = Environment(loader=FileSystemLoader("/workspaces/Epitopedia/epitopedia/viz/templates"))
#         template = env.get_template("index.html")
#         output_from_parsed_template = template.render(data=data["results"].items(), parameters=data["parameters"])
#         handle.write(output_from_parsed_template)


with open("/app/output/" +os.path.basename(os.environ["EPITOPEDIA_DATA_DIR"])) as input_handle:
    data = json.load(input_handle)


app = Flask(
    __name__,
    template_folder="/app/epitopedia/viz/templates",
    static_url_path="/viz/static",
    static_folder="/app/epitopedia/viz/motif_align_viz_js/",
)

@app.get("/")
def main():
    return render_template("index.html", data=data["results"].items(), parameters=data["parameters"])

@app.get("/viz/<hit>")
def viz(hit):
    hit = hit.split(",")
    hit_data = data["results"][hit[0]][int(hit[1])]
    return render_template(
        "viz.html",
        data=hit_data,
        motif_data=zip(hit[0], hit_data["SeqBMM Input Struc Res Nums"], hit_data["EPI_PDB Rep Res Nums"]),
        query_chain=hit_data["EPI_SEQ Input Structure"].split("_")[-1],
        target_chain=hit_data["EPI_PDB Rep PDB"].split("_")[-1],
    )



@app.get("/viz/cif/<id>")
def get_cif(id):
    id = id.removesuffix(".cif").rsplit("_", 1)[0]
    if id.startswith("AF-"):

        doc = cif.read_file(f"{config.AFDB_DIR}/{id}.cif")
        block = doc.sole_block()
        new_doc = cif.Document()
        new_block = new_doc.add_new_block(block.name)
        new_block.add_item(block.find_loop_item("_atom_site.group_PDB"))

        response = make_response(new_doc.as_string())
        response.mimetype = "text/plain"
        return response
    else:
        response = make_response(open(f"{config.PDB_DATABASE_DIR}/{id[1:3].lower()}/{id.lower()}.cif").read())
        response.mimetype = "text/plain"
        return response

@app.get("/viz/aln/<path>")
def get_aln(path):

    response = make_response(open(f"{config.TMALIGN_DIR}/{path}").read())
    response.mimetype = "text/plain"
    return response

@app.get("/figure/rmsd/<hit>")
def get_rmsd_figure(hit):
    hit = hit.split(",")

    return send_file(f"{config.FIGURE_DIR}/rmsd_{hit[0]}_{hit[1]}.png")

@app.get("/figure/episcore/<hit>")
def get_episcore_figure(hit):
    hit = hit.split(",")

    return send_file(f"{config.FIGURE_DIR}/episcore_{hit[0]}_{hit[1]}.png")
   
    # app.run(host="0.0.0.0", ) # debug=True, use_reloader=True
    # Flask.run(host="0.0.0.0",debug=False)