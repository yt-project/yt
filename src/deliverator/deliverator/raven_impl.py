from turbozsi import wsexpose
import turbogears
from generated.deliverator_services_server import DeliveratorServer
import model
import base64

import logging
log = logging.getLogger("deliverator.controllers")

AcceptedKeys = [ # Add keys here
               ]

class DeliveratorImpl(DeliveratorServer):

    @wsexpose()
    def SubmitNewImage(self, req, resp):
        # For now we will just submit it, without question.
        # This is stupid.  Someone, go back in time and slap Previous Matt in
        # the face.
        if not req.APIKey in AcceptedKeys:
            log.debug("APIKey invalid: %s" % (req.APIKey))
            resp.Response = "API Key Invalid"
            return resp
        pf = model.ParameterFile.select( \
             model.ParameterFile.q.GeneratedAt == \
             req.ParameterFile)[0]
        r = model.EnzoRun.select( \
             model.EnzoRun.q.id ==
             req.RunID)[0]
        print "Thingie"
        p = model.Image( \
                parameterfile=pf,
                Field1=str(req.Field1),
                Field2=str(req.Field2),
                Field3=str(req.Field3),
                Axis=str(req.Axis),
                Width=req.Width,
                Unit=str(req.Unit),
                IMG_src=str(req.IMG_src),
                metaData=str(req.MetaData),
                Type=str(req.Type),
                enzorun = r)
        print "Thingie 2!"
        model.hub.commit()
        resp.Response = "Image submitted"
        log.debug("Image submitted with APIKey %s at %s" % (req.APIKey, req.IMG_src))
        return resp

    @wsexpose()
    def QueryExistingRuns(self, req, resp):
        # We only return the first one, which is stupid, but acceptable
        if not req.APIKey in AcceptedKeys:
            log.debug("APIKey invalid: %s" % (req.APIKey))
            resp.RunID = -1
            resp.Response = "API Key Invalid"
            return resp
        r=list(model.EnzoRun.select(model.EnzoRun.q.metaData == req.MetaData))
        if len(r) == 0:
            resp.RunID = -1
            resp.Response = "No run found!"
        else:
            resp.RunID = r[0].id
            resp.Response = "Run found"
        model.hub.commit()
        return resp

    @wsexpose()
    def SubmitNewRun(self, req, resp):
        # For now we will just submit it, without question.
        # This is stupid.  Someone, go back in time and slap Previous Matt in
        # the face.
        if not req.APIKey in AcceptedKeys:
            log.debug("APIKey invalid: %s" % (req.APIKey))
            resp.RunID = -1
            resp.Response = "API Key Invalid"
            return resp
        print "Getting runs"
        r=list(model.EnzoRun.select(model.EnzoRun.q.metaData == str(req.MetaData)))
        print "Got runs"
        if len(r) == 0:
            r = model.EnzoRun(user=str(req.User), metaData=str(req.MetaData))
            resp.RunID = r.id
            resp.Response = "Run submitted"
            model.hub.commit()
        else:
            resp.RunID = r[0].id
            resp.Response = "Run found"
        log.debug("Run submitted with APIKey %s and RunID %s" % (req.APIKey, resp.RunID))
        return resp

    @wsexpose()
    def SubmitNewParameterFile(self, req, resp):
        # For now we will just submit it, without question.
        # This is stupid.  Someone, go back in time and slap Previous Matt in
        # the face.
        if not req.APIKey in AcceptedKeys:
            log.debug("APIKey invalid: %s" % (req.APIKey))
            resp.Response = "API Key Invalid"
            return resp
        print "Querying"
        r=list(model.ParameterFile.select(\
            model.ParameterFile.q.GeneratedAt == req.GeneratedAt))
        print "Queried"
        if len(r) == 0:
            print "Getting run"
            r = model.EnzoRun.select( \
                model.EnzoRun.q.id ==
                req.RunID)[0]
            print "Got run"
            newParam = model.ParameterFile( \
                GeneratedAt=int(req.GeneratedAt), enzorun=r,
                metaData=str(req.MetaData), EnzoHierarchy = base64.b64decode(str(req.PickleObj)), \
                FileName = str(req.FileName) )
            print "Got submission"
            resp.Response = "ParameterFile submitted"
            model.hub.commit()
        else:
            resp.Response = "Parameter file already exists"
        log.debug("Parameter file submitted with APIKey %s and GeneratedAt %s" % (req.APIKey, req.GeneratedAt))
        return resp
