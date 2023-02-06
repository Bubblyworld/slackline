# Has to be run directly with python3 for some flask-related reason.

from flask_restful import Resource, Api
from flask_cors import CORS
from flask import Flask, request
from src.api import Constraints, Rig

app = Flask(__name__)
CORS(app)

class Rig(Resource):
    def post(self):
        constraints = Constraints.from_json(request.json)
        rig = constraints.rig()
        return rig.to_dict()

api = Api(app)
api.add_resource(Rig, '/rig')
if __name__ == '__main__':
    app.run(debug=True)
