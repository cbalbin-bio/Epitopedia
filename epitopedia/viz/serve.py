from flask import Flask, render_template, make_response
from jinja2 import Environment, FileSystemLoader
from rich import print
from werkzeug.wrappers import response
import os

from epitopedia.app import config



def write_html(output_path, data):
    with open(output_path, "w") as handle:

        env = Environment(loader=FileSystemLoader("/workspaces/Epitopedia/epitopedia/viz/templates"))
        template = env.get_template("index.html")
        output_from_parsed_template = template.render(data=data.items())
        handle.write(output_from_parsed_template)


def serve_html(data):
    app = Flask(__name__, template_folder="/workspaces/Epitopedia/epitopedia/viz/templates",static_url_path="/viz/static",static_folder="/workspaces/Epitopedia/epitopedia/viz/motif_align_viz_js/")
    @app.get("/")
    def main():
        return render_template("index.html", data=data.items())


    @app.get("/viz/<hit>")
    def viz(hit):
        hit = hit.split(",")
        hit_data = data[hit[0]][int(hit[1])]
        return render_template("viz.html", data=hit_data, motif_data=zip(hit[0],hit_data["SeqBMM Input Struc Res Nums"],hit_data["EPI_PDB Rep Res Nums"]),query_chain=hit_data["EPI_SEQ Input Structure"].split("_")[-1],target_chain=hit_data["EPI_PDB Rep PDB"].split("_")[-1])

    print("[bold green]View results in browser at http://0.0.0.0:5000[/bold green]")

    @app.get("/viz/cif/<id>")
    def get_cif(id):
        id = id.removesuffix(".cif").rsplit("_",1)[0]
        if id.startswith("AF-"):
            
            response = make_response(open(f"{config.AFDB_DIR}/{id}.cif").read())
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

    app.run(host="0.0.0.0", debug=True)







