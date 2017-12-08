from flask import Flask, request
from flask_restful import Resource, Api
from flask.ext.jsonpify import jsonify
import uuid
import os
import chemkin.nasa
import chemkin.parser
import numpy as np


class Session(Resource):
    def post(self):
        sid = str(uuid.uuid1())
        data = request.json['data']
        folder = "/tmp/chemkin/webserver/{}".format(sid)
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, "data.xml"), "w") as f:
            f.write(data)
        return {'status': 'success', 'id': sid}


class Rates(Resource):
    def post(self, sid):
        folder = "/tmp/chemkin/webserver/{}".format(sid)
        nasa = chemkin.nasa.NASACoeffs()
        # create a data parser class
        data_parser = chemkin.parser.DataParser()
        # parse the data file and return an instance of ReactionData class
        reaction_data = data_parser.parse_file(os.path.join(folder, "data.xml"), nasa)

        conc = [0] * len(reaction_data.species)

        print(request.json)

        for i, sp in enumerate(reaction_data.species):
            conc[i] = float(request.json[sp])

        T = float(request.json['_temp'])

        progress_rates = reaction_data.get_progress_rate(conc, T)  # type: np.ndarray
        reaction_rates = reaction_data.get_reaction_rate(progress_rates)
        ks = reaction_data.get_k(T)
        result = {
            'progress_rates': progress_rates.tolist(),
            'reaction_rates': reaction_rates.tolist(),
            'ks': ks.tolist(),
            'species': reaction_data.species,
        }
        return jsonify(result)


class WebServer:
    def __init__(self, port):
        self.port = port
        self.app = Flask("chemkin web server")
        self.api = Api(self.app)
        self.api.add_resource(Session, '/session')
        self.api.add_resource(Rates, '/rates/<sid>')

    def start(self):
        self.app.run(port=self.port)
