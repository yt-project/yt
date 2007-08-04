import model
from sqlobject import sqlbuilder
import cherrypy
import os.path

def getValsRID():
    if cherrypy.request.params.has_key("rid") == True:
        t = int(cherrypy.request.params['rid'])
    else:
        t = 0
    return t

def getValsCol(column):
    i = sqlbuilder.table.Image
    c = sqlbuilder.table.Image.__getattr__(column)
    if cherrypy.request.params.has_key("rid") == True:
        rid = int(cherrypy.request.params['rid'])
        where = (i.enzorun_ID == rid)
    else: where = None
    a = sqlbuilder.Select(c, where=where, distinct=True)
    a = model.Image._connection.queryAll(sqlbuilder.sqlrepr(a))
    a = map(str, [r[0] for r in a])
    a.sort()
    return zip(a,a)

def getValsParameterFile():
    i = sqlbuilder.table.ParameterFile
    if cherrypy.request.params.has_key("rid") == True:
        rid = int(cherrypy.request.params['rid'])
        a = model.ParameterFile.select(i.enzorun_ID == rid, distinct=True)
    else: a = model.ParameterFile.select()
    k = [(rec.id,os.path.basename(str(rec.FileName))) for rec in a]
    k.sort()
    return k
    

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
