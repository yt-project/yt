import model
from sqlobject import sqlbuilder
import cherrypy

def getValsRID():
    if cherrypy.request.params.has_key("rid") == True:
        t = int(cherrypy.request.params['rid'])
    else:
        t = 0
    return t

def getValsCol(column):
    c = model.Image.q.__getattr__(column)
    query = "SELECT DISTINCT %s from %s" % \
            (sqlbuilder.sqlrepr(c), model.Image.sqlmeta.table)
    if cherrypy.request.params.has_key("rid") == True:
        query += " WHERE Image.enzorun_id = %s" % \
            (int(cherrypy.request.params['rid']))
    vv = []
    vals = model.Image._connection.queryAll(query)
    if len(vals)==1 and vals[0][0] == None:
        return vv
    for val in vals:
        v = str(val[0])
        vv.append((v,v))
    vv.sort()
    return vv

def getValsParameterFile():
    column = "parameterfileID"
    c = model.Image.q.__getattr__(column)
    query = "SELECT DISTINCT %s from %s" % \
            (sqlbuilder.sqlrepr(c), model.Image.sqlmeta.table)
    if cherrypy.request.params.has_key("rid") == True:
        query += " WHERE Image.enzorun_id = %s" % \
            (int(cherrypy.request.params['rid']))
    query += " ORDER BY %s" % (sqlbuilder.sqlrepr(c))
    vv = []
    vals = model.Image._connection.queryAll(query)
    vals=list(model.ParameterFile.select(sqlbuilder.IN(model.ParameterFile.q.id,vals), distinct=True))
    if len(vals)==1:
        return vv
    for val in vals:
        v1 = val.id
        v2 = str(val.FileName).split("/")[-1]
        vv.append((v2,v1))
    vv.sort()
    vv2 = []
    for v1, v2 in vv:
        vv2.append((v2,v1))
    del vv
    vv = vv2
    return vv

def getValsWidth():
    return getValsCol("Width")

def getValsUnit():
    return getValsCol("Unit")

def getValsField1():
    return getValsCol("Field1")

def getValsField2():
    return getValsCol("Field2")

def getValsField3():
    return getValsCol("Field3")

def getValsAxis():
    return getValsCol("Axis")

def getValsType():
    return getValsCol("Type")

def getValsRuns():
    vv = []
    for run in model.EnzoRun.select():
        r = "(%s) %s" % \
            (run.user, run.metaData)
        vv.append((run.id,r))
    return vv
