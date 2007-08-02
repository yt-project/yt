import model
from sqlobject import sqlbuilder
import cherrypy

def getValsRID():
    if cherrypy.request.params.has_key("rid") == True:
        t = int(cherrypy.request.params['rid'])
    else:
        t = 0
    #return [(t,t)]
    print "RETURNING ", t
    return t

def getValsCol(column):
    ##print "GETTING VALS",column
    c = model.Image.q.__getattr__(column)
    query = "SELECT DISTINCT %s from %s" % \
            (sqlbuilder.sqlrepr(c), model.Image.sqlmeta.table)
    ##print cherrypy.request.params.has_key("rid")
    if cherrypy.request.params.has_key("rid") == True:
        # Now we add a 'where' column
        query += " WHERE Image.enzorun_id = %s" % \
            (int(cherrypy.request.params['rid']))
    #print query
    #query += " ORDER BY %s" % (sqlbuilder.sqlrepr(c))
    vv = []
    vals = model.Image._connection.queryAll(query)
    if len(vals)==1 and vals[0][0] == None:
        return vv
    #print vals
    for val in vals:
        v = str(val[0])
        vv.append((v,v))
    #print column, vv
    vv.sort()
    return vv

def getValsDD():
    column = "parameterfileID"
    #print "DD GETTING VALS",column
    c = model.Image.q.__getattr__(column)
    #print "Hello!"
    query = "SELECT DISTINCT %s from %s" % \
            (sqlbuilder.sqlrepr(c), model.Image.sqlmeta.table)
    ##print cherrypy.request.params.has_key("rid")
    #print "QUERY:", query
    if cherrypy.request.params.has_key("rid") == True:
        # Now we add a 'where' column
        query += " WHERE Image.enzorun_id = %s" % \
            (int(cherrypy.request.params['rid']))
    #print "QUERY2:", query
    query += " ORDER BY %s" % (sqlbuilder.sqlrepr(c))
    vv = []
    vals = model.Image._connection.queryAll(query)
    vals=list(model.ParameterFile.select(sqlbuilder.IN(model.ParameterFile.q.id,vals), distinct=True))
    if len(vals)==1:
        return vv
    #print vals
    for val in vals:
        v1 = val.id
        v2 = str(val.FileName).split("/")[-1]
        vv.append((v2,v1))
    #print column, vv
    vv.sort()
    vv2 = []
    for v1, v2 in vv:
        vv2.append((v2,v1))
    del vv
    vv = vv2
    return vv
    #return getValsCol("parameterfileID")

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
