from flask import Flask, render_template
import json

with open("/Users/christianbalbin/bioinformatics/docker/6XR8_A_best_per_source_seq.json") as input_handle:
    data = json.load(input_handle)


app = Flask(__name__)


@app.get("/")
def main():
    return render_template("index.html", data=data.items())
