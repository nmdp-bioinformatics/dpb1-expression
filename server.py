#!env/bin/python
#
# Copyright (c) 2024 NMDP.
#
# This file is part of DPB1 Expression 
# (see https://github.com/nmdp-bioinformatics/dpb1-expression).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
from flask import Flask, jsonify, request
from flask_restx import Api, Resource, fields
from hla_seq_db.hla_seq_db import HlaDB
from pyard import ARD

version = '3500'
ard = ARD(data_dir='tmp/py-ard')
hla_db = HlaDB(db_version=version, loci=['DPB1'], ard=ard, verbose=True)

app = Flask(__name__)
app.config['JSON_SORT_KEYS'] = False
api = Api(app=app)
ns_type = api.namespace('assignment', description='Assign the expression levels of HLA-DPB1 alleles and genotypes.')
ns_match = api.namespace('matching', description='Evaluating DPB1 matching when considering TCE (T-cell epitope) groups and expression levels.')
ns_sort = api.namespace('sorting', description='Evaluating DPB1 sorting when considering TCE (T-cell epitope) groups and expression levels.')

allele_model = api.model('allele',
    {'allele' : fields.String(
        example="DPB1*04:01:01:01/DPB1*04:01:01:02/DPB1*04:01:01:03/DPB1*04:01:01:04/DPB1*04:01:01:05/DPB1*04:01:01:06/DPB1*04:01:01:07/DPB1*04:01:01:08/DPB1*04:01:01:09/DPB1*04:01:01:10/DPB1*04:01:01:11/DPB1*04:01:01:12/DPB1*04:01:01:13/DPB1*04:01:01:14/DPB1*04:01:01:15/DPB1*04:01:01:16/DPB1*04:01:01:17/DPB1*04:01:01:18/DPB1*04:01:01:19/DPB1*04:01:01:20/DPB1*04:01:01:21/DPB1*04:01:01:22/DPB1*04:01:01:23/DPB1*04:01:01:25/DPB1*04:01:01:26/DPB1*04:01:01:24N",
        required=True)
    #  'sire' : fields.String(
    #      example="API", 
    #      description="Self identified race & ethnicity (SIRE). Choices include AFA (African American), API (Asian Pacific Islander), EURO (European), HIS (Hispanic), MENA (Middle Eastern/N African Coast), NAM (Native American), and UNK.",
    #      required=False)
    })
genotype_model = api.model('genotype',
    {'genotype' : fields.String(
        example="DPB1*01:AETTA+DPB1*03:ACMGK",
        required=True),
     'sire' : fields.String(example="CAU", required=False)})

@ns_type.route('/allele/<string:allele>/sire/<sire>')
@ns_type.route('/allele/<string:allele>')
class Dpb1AlleleGet(Resource):
    def get(self, allele : str, sire=None):
        """
        returns expression annotation of HLA-DPB1 allele(s).
        """
        try:
            annotation = hla_db.annotate_allotype(allele).serialize() #, sire=sire)
            return jsonify(annotation)
        except Exception as e:
            return e.__dict__, 500
    
@ns_type.route('/allele')
class Dpb1AllelePost(Resource):
    @api.expect(allele_model)
    def post(self):
        """
        returns expression annotation of HLA-DPB1 allele(s) (example provided; use this if '/' characters present). 
        Exampled provided. Use this if inputting GL (genotype list) strings with '/' characters.
        The 'sire' (self-identified race & ethnicity) parameter is optional.
        """
        try:
            allele = request.json['allele']
            annotation = hla_db.annotate_allotype(allele).serialize()
            return jsonify(annotation)
        except Exception as e:
            return e.__dict__, 500

@ns_type.route("/genotype/<string:genotype>/sire/<sire>")
@ns_type.route("/genotype/<string:genotype>")
class DPB1GenotypeGet(Resource):
    def get(self, genotype : str, sire : str = None):
        """
        returns expression annotation of HLA-DPB1 genotype(s).
        """
        try:
            annotation = hla_db.annotate_genotype(genotype).serialize()
            return jsonify(annotation)
        except Exception as e:
            return e.__dict__, 500

@ns_type.route("/genotype")
class DPB1GenotypePost(Resource):
    @api.expect(genotype_model)
    def post(self):
        """
        returns expression annotation of HLA-DPB1 genotype(s)  (example provided; use this if '/' characters present).
        Exampled provided. Use this if inputting GL (genotype list) strings with '/' characters.
        """
        try:
            genotype = request.json['genotype']
            annotation = hla_db.annotate_genotype(genotype).serialize()
            return jsonify(annotation)
        except Exception as e:
            return e.__dict__, 500

@ns_type.route('/map')
class Dpb1AllelePost(Resource):

    def get(self):
        """
        returns map of high-resolution (2-field) HLA-DPB1 alleles to predicted expression levels.
        """
        try:
            dp_map = hla_db.get_map('DPB1', 'expression',
               group_col='allele_two_field',
               unk_val='unknown')
            return jsonify(dp_map)
        except Exception as e:
            return e.__dict__, 500

@ns_match.route("/genotypes")
class Match(Resource):
    
    model = api.model('HLA-DPB1 Genotypes Input', 
		  {'dpb1_genotype_patient': fields.String(required = True, 
					 description="HLA-DPB1 genotype of the patient", 
					 help="HLA-DPB1 genotype format: DPB1*XX:XX+DPB1*YY:XX",
                     example="DPB1*04:01+DPB1*40:01"),
            'dpb1_genotype_donors': fields.String(required = True, 
					 description="HLA-DPB1 genotype of the donor", 
					 help="HLA-DPB1 genotype format: DPB1*XX:XX+DPB1*YY:XX",
                     example="DPB1*40:01+DPB1*40:01")})

    @api.expect(model)
    def post(self):
        """
        Returns annotation on DPB1 matching via the TCE (T-Cell Epitope) and expression model  (example provided).
        """
        try:
            geno_recip = request.json['dpb1_genotype_patient']
            geno_donor = request.json['dpb1_genotype_donors']
            sorting_methods = ['expression', 'tce', 'directionality']
            if 'sorting_methods' in request.json:
                sorting_methods = request.json['sorting_methods']
            match_list = hla_db.annotate_match_list(geno_recip, geno_donor)
            match_list.sort(sorting_methods)
            return jsonify(match_list.serialize())
        except Exception as e:
            print(e)
            return e.__dict__, 500

@ns_sort.route("/genotypes")
class Sort(Resource):
    model_pairs = api.model('Pairs', 
            {'dpb1_genotype_patient': fields.String(example="DPB1*02:01+DPB1*03:01"),
             'dpb1_genotype_donors' : fields.List(fields.String,  
                example=['DPB1*03:01+DPB1*04:01',
                            'DPB1*03:01+DPB1*03:ACMGK',
                            'DPB1*03:01+DPB1*17:01', 
                            'DPB1*02:01+DPB1*14:01',
                            'DPB1*02:01+DPB1*09:01',
                            'DPB1*02:01+DPB1*04:01',
                            'DPB1*02:01+DPB1*02:01'
             ])
            })
    @api.expect(model_pairs)
    def post(self):
        """
        Returns annotation on DPB1 sorting via the TCE (T-Cell Epitope) and expression model (example provided).
        """
        try:
            geno_recip = request.json['dpb1_genotype_patient']
            geno_donor = request.json['dpb1_genotype_donors']
            sorting_methods = ['expression', 'tce', 'directionality']
            if 'sorting_methods' in request.json:
                sorting_methods = request.json['sorting_methods']
            match_list = hla_db.annotate_match_list(geno_recip, geno_donor)
            match_list.sort(sorting_methods)
            return jsonify(match_list.serialize())
        except Exception as e:
            print(e)
            return e.__dict__, 500

@ns_sort.route("/ids")
class Sort(Resource):
    subject = api.model('Subject', {
        'id' : fields.String,
        'genotype' : fields.String,
        'sire' : fields.String
    })
    model_pair_ids = api.model('Pairs', 
            {'dpb1_genotype_patient': fields.String(example="DPB1*02:01+DPB1*03:01"),
             'dpb1_genotype_donors' : fields.List(fields.Nested(subject), 
                    example=[{'id' : 'Donor #1', 'genotype' : 'DPB1*02:01+DPB1*04:01', 'sire' : 'EURO'},
                             {'id' : 'Donor #2', 'genotype' : 'DPB1*03:01+DPB1*03:ACMGK', 'sire' : 'API'},
                             {'id' : 'Donor #3', 'genotype' : 'DPB1*02:01+DPB1*09:01', 'sire' : 'HIS'},
                             {'id' : 'Donor #4', 'genotype' : 'DPB1*02:01+DPB1*02:01', 'sire' : 'UNK'},
                             {'id' : 'Donor #5', 'genotype' : 'DPB1*03:01+DPB1*04:01', 'sire' : 'MENAM'},
                             {'id' : 'Donor #6', 'genotype' : 'DPB1*02:01+DPB1*14:01', 'sire' : 'NAM'},
                             {'id' : 'Donor #7', 'genotype' : 'DPB1*03:01+DPB1*17:01'}])
            })
    @api.expect(model_pair_ids)
    def post(self):
        """
        Returns annotation on DPB1 sorting via the TCE (T-Cell Epitope) and expression model (example provided).
        """
        try:
            geno_recip = request.json['dpb1_genotype_patient']
            geno_donor = request.json['dpb1_genotype_donors']
            sorting_methods = None
            if 'sorting_methods' in request.json:
                sorting_methods = request.json['sorting_methods']
            annotation = hla_db.annotate_match_list(geno_recip, geno_donor)
            # print(annotation, sorting_methods)
            if sorting_methods:
                annotation.sort(sorting_methods)
            else:
                annotation.sort()
            return jsonify(annotation.serialize())
        except Exception as e:
            print(e)
            return e.__dict__, 500

@app.after_request
def after_request(response):
    response.headers.add('Access-Control-Allow-Origin', '*')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
    response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE')
    response.headers['Strict-Transport-Security'] = 'max-age=31536000; includeSubDomains'
    response.headers['Content-Security-Policy'] = "default-src 'self' 'unsafe-inline'"
    response.headers['X-Content-Type-Options'] = 'nosniff'
    response.headers['X-XSS-Protection'] = '1; mode=block'

    return response

if __name__ == '__main__':
    app.run("0.0.0.0", "5010", debug=True)