from cherrypy import request, response
# from deliverator import json
# import logging
# log = logging.getLogger("deliverator.controllers")

from turbogears import controllers, expose, validate, validators, identity, redirect, widgets
from turbogears.toolbox.catwalk import CatWalk
import model, types
from sqlobject.sqlbuilder import Select
from sqlobject import sqlbuilder, func
from raven_impl import DeliveratorImpl
from getoptionvals import *

import cherrypy, string

import logging
log = logging.getLogger("deliverator.controllers")

convDict = {
    "Width": (float, int, func.ROUND)
    }

def NoneFunc(data):
    return data

class MultipleSelectFieldSelectAll(widgets.MultipleSelectField):
    template = '''
    <table>
    <tr>
    <th rowspan="2">
        <select xmlns:py="http://purl.org/kid/ns#"
            multiple="multiple"
            size="${size}"
            name="${name}"
            class="${field_class}"
            id="${field_id}"
            py:attrs="attrs"
        >   
            <optgroup py:for="group, options in grouped_options"
                label="${group}"
                py:strip="not group"
            >   
                <option py:for="value, desc, attrs in options"
                    value="${value}"
                    py:attrs="attrs"
                    py:content="desc"
                />
            </optgroup>
        </select>
    </th>
    <td>
        <input type="button" value="All"
onClick="javascript:SelectAllList(this.form.${name});update_count();"/><br/> <input type="button" value="None" onClick="javascript:DeSelectAllList(this.form.${name});update_count();"/>
    </td></tr>
    <tr><td>
    &nbsp;
    </td></tr>
    </table>
'''

class EnzoRunIDHiddenField(widgets.HiddenField):
    template = '''
    <?python
    import cherrypy as http
    ?>
    <input xmlns:py="http://purl.org/kid/ns#"  
        type="hidden"  
        name="${name}"  
        class="${field_class}"  
        id="${field_id}"  
        value="${str(attrs)}"
    />  
'''

class ImageSelectionForm(widgets.WidgetsList):
    parameterfileID = MultipleSelectFieldSelectAll(validator=validators.Int, \
        options=getValsDD, defaults=None, attrs={'onchange':'js:update_count();'})
    Type = MultipleSelectFieldSelectAll(validator=validators.String, \
        options=getValsType, defaults=None, attrs={'onchange':'js:update_count();'})
    Width = MultipleSelectFieldSelectAll(validator=validators.Int, \
        options=getValsWidth, default=None, attrs={'onchange':'js:update_count();'})
    Unit = MultipleSelectFieldSelectAll(validator=validators.String, \
        options=getValsUnit, default=None, attrs={'onchange':'js:update_count();'})
    Field1 = MultipleSelectFieldSelectAll(validator=validators.String, \
        options=getValsField1, defaults=None, attrs={'onchange':'js:update_count();'})
    Field2 = MultipleSelectFieldSelectAll(validator=validators.String, \
        options=getValsField2, defaults=None, attrs={'onchange':'js:update_count();'})
    Field3 = MultipleSelectFieldSelectAll(validator=validators.String, \
        options=getValsField3, defaults=None, attrs={'onchange':'js:update_count();'})
    Axis = MultipleSelectFieldSelectAll(validator=validators.String, \
        options=getValsAxis, defaults=None, attrs={'onchange':'js:update_count();'})

class RunSelectionForm(widgets.WidgetsList):
    rid = widgets.SingleSelectField(validator=validators.Int, \
        options = getValsRuns, default=None, name="rid")

class Root(controllers.RootController):
    @expose(template="deliverator.templates.welcome")
    # @identity.require(identity.in_group("admin"))
    def index(self):
        raise redirect("/Deliverator/selectRun")
        import time
        log.debug("Happy TurboGears Controller Responding For Duty")
        return dict(now=time.ctime())


    @expose(template="deliverator.templates.login")
    def login(self, forward_url=None, previous_url=None, *args, **kw):

        if not identity.current.anonymous \
            and identity.was_login_attempted() \
            and not identity.get_identity_errors():
            #raise redirect(forward_url)
            raise redirect("/Deliverator/")

        forward_url=None
        previous_url= request.path

        if identity.was_login_attempted():
            msg=_("The credentials you supplied were not correct or "
                   "did not grant access to this resource.")
        elif identity.get_identity_errors():
            msg=_("You must provide your credentials before accessing "
                   "this resource.")
        else:
            msg=_("Please log in.")
            forward_url= request.headers.get("Referer", "/Deliverator/")
            
        response.status=403
        return dict(message=msg, previous_url=previous_url, logging_in=True,
                    original_parameters=request.params,
                    forward_url=forward_url)

    @expose()
    def logout(self):
        identity.current.logout()
        raise redirect("/Deliverator/")

    @expose(template="deliverator.templates.welcome")
    def index(self):
        import time
        # log.debug("Happy TurboGears Controller Responding For Duty")
        return dict(now=time.ctime())

    #@expose()
    #catwalk = CatWalk(model)

    @expose(template="deliverator.templates.selection")
    #@identity.require(identity.in_group('users'))
    def selectImages(self, *a, **kw):
        if kw.has_key("rid"):
            rid=kw["rid"]
        else:
            raise redirect("/Deliverator/selectRun")
        enzorunID = widgets.HiddenField(attrs={'value':rid}, name="enzorunID", defaults=None)
        select_form = widgets.TableForm(fields = ImageSelectionForm(), hidden_fields=[enzorunID])
        vals = {'form_enzorunID':rid}
        return dict(form=select_form, action='gallery', rid=rid, submit_text='Query', title='New Comment')

    @expose(template="deliverator.templates.welcome")
    def moron(self, filename="testgallery.pkl"):
        import pickle, time
        a = pickle.load(open(filename))
        for image in a:
            j = {}
            j["IMG_src"] = image["Filename"]
            j["Width"] = image["Width"]
            j["Unit"] = image["Unit"]
            j["Axis"] = image["Axis"]
            j["Field1"] = image["Field"]
            j["Field2"] = None
            j["Field3"] = None
            j["metaData"] = None
            i = model.Images(**j)
        return dict(now=time.ctime())

    #@validate(ImageQuery)
    def generateQuery(self, data):
        # We're going to be querying and grabbing
        k = data.keys()
        k.sort()
        s = []
        if len(k) == 0:
            return None
        i = 0
        for key in k:
            v = data[key]
            i += 1
            c = model.Image.q.__getattr__(key)
            #print key, v, type(v)
            if not isinstance(v, types.ListType):
                v = [v]
            if convDict.has_key(key):
                v = map(convDict[key][0],v)
                v = map(convDict[key][1],v)
                c = convDict[key][2](c)
            else:
                v = map(str,v)
            options = ",".join(map(repr,v))
            s.append("%s IN (%s)" % (sqlbuilder.sqlrepr(c), options))
        ss="(" + ") AND (".join(s) + ")"
        #print ss
        if i == 0:
            return None
        return ss

    @expose(template="deliverator.templates.gallery")
    def gallery(self, **data):
        # We assume the gallery is only generated for a single Run
        # Not a bad assumption...
        ss = self.generateQuery(data)
        sr = model.Image.select(ss)
        ll = list(sr)
        show_delete = False
        if len(ll) > 0:
            # User ...
            s = ll[0].enzorun.user
            if identity.current.user:
                show_delete = (s == identity.current.user.user_name)
        return dict(images=ll, num=len(ll), show_delete = show_delete)

    @expose(allow_json=True)
    #@validate(validators=dict(Width=validators.Int()))
    def countImages(self, **data):
        dd = {}
        for key in data.keys():
            # Now we strip out the first element
            dd[key] = []
            #print data[key].split(",")[1:]
            for l in data[key].split(",")[1:]:
                if len(l) > 0:
                    dd[key].append(l)
            if dd[key] == []:
                del dd[key]
            #dd[key] = map(string.strip,data[key].split(",")[1:])
        ss = self.generateQuery(dd)
        sr = model.Image.select(ss)
        cc=sr.count()
        #print "COUNT", cc
        return "%s" % (cc)

    @expose(template="deliverator.templates.welcome")
    def deleteImages(self, **data):
        if not identity.current.user or not data.has_key("todelete"):
            raise redirect("/Deliverator/selectRun")
        images_to_verify = map(int,data["todelete"])
        images_to_delete = []
        for image in model.Image.select(sqlbuilder.IN(model.Image.q.id, images_to_verify)):
            if image.enzorun.user == identity.current.user.user_name:
                model.Image.delete(image.id)
        #model.Image.delete(sqlbuilder.IN(model.Image.q.id, images_to_delete))
        print request.path
        raise redirect("/Deliverator/selectRun")

    @expose(template="deliverator.templates.welcome")
    def deleteRun(self, **data):
        if not identity.current.user or not data.has_key("todelete"):
            raise redirect("/Deliverator/selectRun")
        rs_to_verify = data["todelete"]
        if not isinstance(rs_to_verify, types.ListType):
            rs_to_verify = [rs_to_verify]
        rs_to_verify = map(int, rs_to_verify)
        rs_to_delete = []
        print "DELETION!", rs_to_verify
        for run in model.EnzoRun.select(sqlbuilder.IN(model.EnzoRun.q.id, rs_to_verify)):
            if run.user == identity.current.user.user_name:
                # Now we delete all the images first, then the pfile
                for pf in run.ParameterFiles:
                    for image in pf.Images:
                        print "Deleting %s" % (image.IMG_src)
                        model.Image.delete(image.id)
                    print "Deleting %s" % (pf.GeneratedAt)
                    model.ParameterFile.delete(pf.id)
                print "Deleting %s" % (run.id)
                model.EnzoRun.delete(run.id)
        #print request.path
        raise redirect("/Deliverator/selectRun")

    @expose(template="deliverator.templates.welcome")
    def deleteParamFiles(self, **data):
        if not identity.current.user or not data.has_key("todelete"):
            raise redirect("/Deliverator/selectRun")
        pfs_to_verify = data["todelete"]
        if not isinstance(pfs_to_verify, types.ListType):
            pfs_to_verify = [pfs_to_verify]
        pfs_to_verify = map(int, pfs_to_verify)
        pfs_to_delete = []
        print "DELETION!", pfs_to_verify
        for pf in model.ParameterFile.select(sqlbuilder.IN(model.ParameterFile.q.GeneratedAt, pfs_to_verify)):
            if pf.enzorun.user == identity.current.user.user_name:
                # Now we delete all the images first, then the pfile
                for image in pf.Images:
                    #print "Deleting %s" % (image.IMG_src)
                    model.Image.delete(image.id)
                #print "Deleting %s" % (pf.GeneratedAt)
                model.ParameterFile.delete(pf.id)
        #print request.path
        raise redirect("/Deliverator/selectRun")

    @expose(template="deliverator.templates.paraminfo")
    def paraminfo(self, **data):
        if not data.has_key("id"):
            raise redirect("/Deliverator/")
        ids = data["id"]
        if not isinstance(ids, types.ListType):
            ids = [ids]
        try:
            ids = map(int, ids)
        except:
            raise redirect("/Deliverator/")
        if identity.current.user:
            myuser = identity.current.user.user_name
        else:
            myuser = None
        pf = list(model.ParameterFile.select(sqlbuilder.IN(model.ParameterFile.q.GeneratedAt, ids)))
        return dict(pfs=pf, myuser=myuser)

    @expose(template="deliverator.templates.listparams")
    def listParams(self, **kw):
        if kw.has_key("rid"):
            rid=int(kw["rid"])
        else:
            raise redirect("/Deliverator/selectRun")
        run = model.EnzoRun.select(model.EnzoRun.q.id == rid)[0]
        show_delete = False
        if identity.current.user:
            show_delete = (run.user == identity.current.user.user_name)
        return dict(run=run, show_delete=show_delete)

    @expose(template="deliverator.templates.selectrun")
    def selectRun(self, **data):
        # We construct a list of Enzo Runs
        r = model.EnzoRun.select()
        select_form = widgets.TableForm(fields = RunSelectionForm())
        return dict(form=select_form, action='selectImages', submit_text='Query', title='TITLE')

    #@expose(template="deliverator.templates.runinfo")
    #def runInfo(self, **data):
        #return ""

    RavenMethods = DeliveratorImpl()
