from cherrypy import request, response
# from deliverator import json
# import logging
# log = logging.getLogger("deliverator.controllers")

from turbogears import controllers, expose, validate, validators, identity, redirect, widgets
from turbogears.toolbox.catwalk import CatWalk
import model, types
from sqlobject.sqlbuilder import Select
from sqlobject import sqlbuilder, func
from deliverator_soap import DeliveratorImpl
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
                    py:content="desc"
                    py:attrs="attrs"
                />
            </optgroup>
        </select>
    </th>
    <td>
        <input type="button" value="All" style="width:50px;"
onClick="javascript:SelectAllList(this.form.${name});update_count();"/><br/>
        <input type="button" value="None" style="width:50px;"
onClick="javascript:DeSelectAllList(this.form.${name});update_count();"/>
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

# There may be an easier way to do this.
# This is the cleanest way I could think of off-hand.
def getSelForm():
    def getMSFSA(name, val, f):
        j = MultipleSelectFieldSelectAll(validator=val, \
            options=f, defaults=None, \
            attrs={'onchange':'js:update_count();'},\
            field_class="msfsa")
        return j
    class ImageSelectionFormBase(widgets.WidgetsList):
        ParameterFile = getMSFSA("ParameterFile", validators.Int, getValsParameterFile)
        Type = getMSFSA("Type", validators.String, getValsType)
        Width = getMSFSA("Width", validators.Int, getValsWidth)
        Unit = getMSFSA("Unit", validators.String, getValsUnit)
        Field1 = getMSFSA("Field1", validators.String, getValsField1)
        Field2 = getMSFSA("Field2", validators.String, getValsField2)
        Field3 = getMSFSA("Field3", validators.String, getValsField3)
        Axis = getMSFSA("Axis", validators.String, getValsAxis)
    return ImageSelectionFormBase

ImageSelectionForm = getSelForm()

class RunSelectionForm(widgets.WidgetsList):
    rid = widgets.SingleSelectField(validator=validators.Int, \
        options = getValsRuns, default=None, name="rid")


class SelectImagesController:
    @expose(template="deliverator.templates.selection")
    def default(self, *a, **kw):
        if kw.has_key("rid"):
            raise redirect("/Deliverator/selectImages/%s" % (kw['rid']))
        else:
            rid=int(a[0])
            print "GOT RID:", rid
        enzorunID = widgets.HiddenField(attrs={'value':rid}, name="enzorunID", defaults=None)
        select_form = widgets.TableForm(fields = ImageSelectionForm(), hidden_fields=[enzorunID])
        vals = {'form_enzorunID':rid}
        return dict(form=select_form, action='/Deliverator/gallery', rid=rid, submit_text='Query', title='New Comment')

class Root(controllers.RootController):

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
        return dict()

    selectImages = SelectImagesController()

    def generateQuery(self, data):
        def ensureList(v):
            if not isinstance(v,types.ListType):
                v = [v]
            return v
        queries = []
        i = sqlbuilder.table.Image
        if data.has_key("ParameterFile"):
            pids = ensureList(data.pop("ParameterFile"))
            queries.append(sqlbuilder.IN(i.parameterfile_ID,
                              map(int, pids)))
        if data.has_key("parameterfileID"):
            pids = ensureList(data.pop("parameterfileID"))
            queries.append(sqlbuilder.IN(i.parameterfile_ID,
                              map(int, pids)))
        if data.has_key("enzorunID"):
            pids = ensureList(data.pop("enzorunID"))
            queries.append(sqlbuilder.IN(i.enzorun_ID,
                              map(int, pids)))
        if data.has_key("Width"):
            widths = ensureList(data.pop("Width"))
            queries.append(sqlbuilder.IN(func.ROUND(i.Width),
                              map(int, map(float, widths))))
        for k, v in data.items():
            v = ensureList(v)
            queries.append(sqlbuilder.IN(i.__getattr__(k),map(str,v)))
        k = sqlbuilder.AND(*queries)
        return k

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
            for l in data[key].split(",")[1:]:
                if len(l) > 0:
                    dd[key].append(l)
            if dd[key] == []:
                del dd[key]
        ss = self.generateQuery(dd)
        sr = model.Image.select(ss)
        cc=sr.count()
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
        for run in model.EnzoRun.select(sqlbuilder.IN(model.EnzoRun.q.id, rs_to_verify)):
            if run.user == identity.current.user.user_name:
                # Now we delete all the images first, then the pfile
                for pf in run.ParameterFiles:
                    for image in pf.Images:
                        model.Image.delete(image.id)
                    model.ParameterFile.delete(pf.id)
                model.EnzoRun.delete(run.id)
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
        for pf in model.ParameterFile.select(sqlbuilder.IN(model.ParameterFile.q.GeneratedAt, pfs_to_verify)):
            if pf.enzorun.user == identity.current.user.user_name:
                # Now we delete all the images first, then the pfile
                for image in pf.Images:
                    model.Image.delete(image.id)
                model.ParameterFile.delete(pf.id)
        raise redirect("/Deliverator/selectRun")

    @expose(template="deliverator.templates.paraminfo")
    def paraminfo(self, **data):
        if not data.has_key("id"): raise redirect("/Deliverator/")
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

    @expose(template="deliverator.templates.viewimage")
    def viewimage(self, **data):
        if not data.has_key("id"): raise redirect("/Deliverator/")
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
        ims = list(model.Image.select(sqlbuilder.IN(model.Image.q.id, ids)))[0]
        return dict(Image=ims, myuser=myuser)

    

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
        return dict(form=select_form, action='chooseRun', submit_text='Query', title='TITLE')

    def chooseRun(self, *a, **kw):
        if not kw.has_key('rid'):
            raise redirect("/Deliverator/selectRun")
        rid = kw['rid']
        raise redirect("/Deliverator/selectImages/%s" % (rid))

    @expose(template="deliverator.templates.ajaxselectrun", allow_json=True)
    def ajaxselectrun(self, **data):
        pass

    DeliveratorMethods = DeliveratorImpl()
