# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.

import json

from flask import Flask, render_template

with open("/Users/christianbalbin/bioinformatics/docker/6XR8_A_best_per_source_seq.json") as input_handle:
    data = json.load(input_handle)


app = Flask(__name__)


def saveHTML(path):
    with open(path, "w") as handle:
        handle.write(render_template("index.html", data=data.items()))


@app.get("/")
def main():
    return render_template("index.html", data=data.items())
