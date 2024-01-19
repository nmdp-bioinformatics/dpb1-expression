from hla_seq_db.hla_seq_db import HlaDB

def before_all(context):
    context.hla_db = HlaDB(db_version='3470', verbose=True)