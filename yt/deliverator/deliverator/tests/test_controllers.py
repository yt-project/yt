import turbogears
from nose import with_setup
from turbogears import testutil
from deliverator.controllers import Root
import cherrypy

def teardown_func():
    """Tests for apps using identity need to stop CP/TG after each test to
    stop the VisitManager thread. See http://trac.turbogears.org/turbogears/ticket/1217
    for details.
    """
    turbogears.startup.stopTurboGears()

cherrypy.root = Root()

def test_method():
    "the index method should return a string called now"
    import types
    result = testutil.call(cherrypy.root.index)
    assert type(result["now"]) == types.StringType
test_method = with_setup(teardown=teardown_func)(test_method)

def test_indextitle():
    "The indexpage should have the right title"
    testutil.createRequest("/")
    assert "<TITLE>Welcome to TurboGears</TITLE>" in cherrypy.response.body[0]
test_indextitle = with_setup(teardown=teardown_func)(test_indextitle)

def test_logintitle():
    "login page should have the right title"
    testutil.createRequest("/login")
    assert "<TITLE>Login</TITLE>" in cherrypy.response.body[0]
test_logintitle = with_setup(teardown=teardown_func)(test_logintitle)
