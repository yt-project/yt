"""
This is a modified version of the post-review script from the package RBTools,
available on PyPI at http://pypi.python.org/pypi/RBTools/ .  The authors of the
file recommended not distributing a modified version, but we're going to try to
keep in sync anyway.  The local modifications were implemented by Matthew Turk
<matthewturk@gmail.com> and include:

    * Removal of some of the version control systems, leaving only Mercurial
      and SVN
    * Some pre-loading of information
    * Some auto-prompting of information (takes place in yt/commands.py)
    * Made options a non-global

The original AUTHORS file contained:

Lead Developers:

	* Christian Hammond
	* David Trowbridge


Contributors:

	* Chris Clark
	* Eric Huss
	* Jeremy Bettis
	* Lepton Wu
	* Luke Lu
	* Paul Scott
	* Raghu Kaippully
	* Stacey Sheldon
	* Steven Russell

and the original license is:

Copyright (c) 2007  Christian Hammond
Copyright (c) 2007  David Trowbridge

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
#!/usr/bin/env python
import cookielib
import difflib
import getpass
import marshal
import mimetools
import ntpath
import os
import re
import socket
import stat
import subprocess
import sys
import tempfile
import urllib
import urllib2
from optparse import OptionParser
from tempfile import mkstemp
from urlparse import urljoin, urlparse
from yt.logger import ytLogger as mylog

try:
    import json
except ImportError:
    import simplejson as json

# This specific import is necessary to handle the paths for
# cygwin enabled machines.
if (sys.platform.startswith('win')
    or sys.platform.startswith('cygwin')):
    import ntpath as cpath
else:
    import posixpath as cpath

###
# Default configuration -- user-settable variables follow.
###

# The following settings usually aren't needed, but if your Review
# Board crew has specific preferences and doesn't want to express
# them with command line switches, set them here and you're done.
# In particular, setting the REVIEWBOARD_URL variable will allow
# you to make it easy for people to submit reviews regardless of
# their SCM setup.
#
# Note that in order for this script to work with a reviewboard site
# that uses local paths to access a repository, the 'Mirror path'
# in the repository setup page must be set to the remote URL of the
# repository.

#
# Reviewboard URL.
#
# Set this if you wish to hard-code a default server to always use.
# It's generally recommended to set this using your SCM repository
# (for those that support it -- currently only SVN, Git, and Perforce).
#
# For example, on SVN:
#   $ svn propset reviewboard:url http://reviewboard.example.com .
#
# Or with Git:
#   $ git config reviewboard.url http://reviewboard.example.com
#
# On Perforce servers version 2008.1 and above:
#   $ p4 counter reviewboard.url http://reviewboard.example.com
#
# Older Perforce servers only allow numerical counters, so embedding
# the url in the counter name is also supported:
#   $ p4 counter reviewboard.url.http:\|\|reviewboard.example.com 1
#
# Note that slashes are not allowed in Perforce counter names, so replace them
# with pipe characters (they are a safe substitute as they are not used
# unencoded in URLs). You may need to escape them when issuing the p4 counter
# command as above.
#
# If this is not possible or desired, setting the value here will let
# you get started quickly.
#
# For all other repositories, a .reviewboardrc file present at the top of
# the checkout will also work. For example:
#
#   $ cat .reviewboardrc
#   REVIEWBOARD_URL = "http://reviewboard.example.com"
#
REVIEWBOARD_URL = "http://review.enzotools.org/"

###
# End user-settable variables.
###

VERSION = "0.8"

user_config = None
tempfiles = []
options = None

class APIError(Exception):
    pass


class RepositoryInfo:
    """
    A representation of a source code repository.
    """
    def __init__(self, path=None, base_path=None, supports_changesets=False,
                 supports_parent_diffs=False):
        self.path = path
        self.base_path = base_path
        self.supports_changesets = supports_changesets
        self.supports_parent_diffs = supports_parent_diffs
        mylog.debug("repository info: %s" % self)

    def __str__(self):
        return "Path: %s, Base path: %s, Supports changesets: %s" % \
            (self.path, self.base_path, self.supports_changesets)

    def set_base_path(self, base_path):
        if not base_path.startswith('/'):
            base_path = '/' + base_path
        mylog.debug("changing repository info base_path from %s to %s" % \
              (self.base_path, base_path))
        self.base_path = base_path

    def find_server_repository_info(self, server):
        """
        Try to find the repository from the list of repositories on the server.
        For Subversion, this could be a repository with a different URL. For
        all other clients, this is a noop.
        """
        return self


class SvnRepositoryInfo(RepositoryInfo):
    """
    A representation of a SVN source code repository. This version knows how to
    find a matching repository on the server even if the URLs differ.
    """
    def __init__(self, path, base_path, uuid):
        RepositoryInfo.__init__(self, path, base_path)
        self.uuid = uuid

    def find_server_repository_info(self, server):
        """
        The point of this function is to find a repository on the server that
        matches self, even if the paths aren't the same. (For example, if self
        uses an 'http' path, but the server uses a 'file' path for the same
        repository.) It does this by comparing repository UUIDs. If the
        repositories use the same path, you'll get back self, otherwise you'll
        get a different SvnRepositoryInfo object (with a different path).
        """
        repositories = server.get_repositories()

        for repository in repositories:
            if repository['tool'] != 'Subversion':
                continue

            info = self._get_repository_info(server, repository)

            if not info or self.uuid != info['uuid']:
                continue

            repos_base_path = info['url'][len(info['root_url']):]
            relpath = self._get_relative_path(self.base_path, repos_base_path)
            if relpath:
                return SvnRepositoryInfo(info['url'], relpath, self.uuid)

        # We didn't find a matching repository on the server. We'll just return
        # self and hope for the best.
        return self

    def _get_repository_info(self, server, repository):
        try:
            return server.get_repository_info(repository['id'])
        except APIError, e:
            # If the server couldn't fetch the repository info, it will return
            # code 210. Ignore those.
            # Other more serious errors should still be raised, though.
            rsp = e.args[0]
            if rsp['err']['code'] == 210:
                return None

            raise e

    def _get_relative_path(self, path, root):
        pathdirs = self._split_on_slash(path)
        rootdirs = self._split_on_slash(root)

        # root is empty, so anything relative to that is itself
        if len(rootdirs) == 0:
            return path

        # If one of the directories doesn't match, then path is not relative
        # to root.
        if rootdirs != pathdirs:
            return None

        # All the directories matched, so the relative path is whatever
        # directories are left over. The base_path can't be empty, though, so
        # if the paths are the same, return '/'
        if len(pathdirs) == len(rootdirs):
            return '/'
        else:
            return '/'.join(pathdirs[len(rootdirs):])

    def _split_on_slash(self, path):
        # Split on slashes, but ignore multiple slashes and throw away any
        # trailing slashes.
        split = re.split('/*', path)
        if split[-1] == '':
            split = split[0:-1]
        return split


class ReviewBoardHTTPPasswordMgr(urllib2.HTTPPasswordMgr):
    """
    Adds HTTP authentication support for URLs.

    Python 2.4's password manager has a bug in http authentication when the
    target server uses a non-standard port.  This works around that bug on
    Python 2.4 installs. This also allows post-review to prompt for passwords
    in a consistent way.

    See: http://bugs.python.org/issue974757
    """
    def __init__(self, reviewboard_url):
        self.passwd  = {}
        self.rb_url  = reviewboard_url
        self.rb_user = None
        self.rb_pass = None

    def find_user_password(self, realm, uri):
        if uri.startswith(self.rb_url):
            if self.rb_user is None or self.rb_pass is None:
                print "==> HTTP Authentication Required"
                print 'Enter username and password for "%s" at %s' % \
                    (realm, urlparse(uri)[1])
                self.rb_user = raw_input('Username: ')
                self.rb_pass = getpass.getpass('Password: ')

            return self.rb_user, self.rb_pass
        else:
            # If this is an auth request for some other domain (since HTTP
            # handlers are global), fall back to standard password management.
            return urllib2.HTTPPasswordMgr.find_user_password(self, realm, uri)


class ReviewBoardServer(object):
    """
    An instance of a Review Board server.
    """
    def __init__(self, url, info, cookie_file, options):
        self.options = options
        self.url = url
        if self.url[-1] != '/':
            self.url += '/'
        self._info = info
        self._server_info = None
        self.cookie_file = cookie_file
        self.cookie_jar  = cookielib.MozillaCookieJar(self.cookie_file)

        # Set up the HTTP libraries to support all of the features we need.
        cookie_handler = urllib2.HTTPCookieProcessor(self.cookie_jar)
        password_mgr   = ReviewBoardHTTPPasswordMgr(self.url)
        auth_handler   = urllib2.HTTPBasicAuthHandler(password_mgr)

        opener = urllib2.build_opener(cookie_handler, auth_handler)
        opener.addheaders = [('User-agent', 'post-review/' + VERSION)]
        urllib2.install_opener(opener)

    def login(self, force=False):
        """
        Logs in to a Review Board server, prompting the user for login
        information if needed.
        """
        if not force and self.has_valid_cookie():
            return

        print "==> Review Board Login Required"
        print "Enter username and password for Review Board at %s" % self.url
        print "(No username?  Get one here: %s/account/register )" % (self.url)
        if self.options.username:
            username = self.options.username
        elif self.options.submit_as:
            username = self.options.submit_as
        else:
            username = raw_input('Username: ')

        if not self.options.password:
            password = getpass.getpass('Password: ')
        else:
            password = self.options.password

        mylog.debug('Logging in with username "%s"' % username)
        try:
            self.api_post('api/json/accounts/login/', {
                'username': username,
                'password': password,
            })
        except APIError, e:
            rsp, = e.args

            die("Unable to log in: %s (%s)" % (rsp["err"]["msg"],
                                               rsp["err"]["code"]))

        mylog.debug("Logged in.")

    def has_valid_cookie(self):
        """
        Load the user's cookie file and see if they have a valid
        'rbsessionid' cookie for the current Review Board server.  Returns
        true if so and false otherwise.
        """
        try:
            parsed_url = urlparse(self.url)
            host = parsed_url[1]
            path = parsed_url[2] or '/'

            # Cookie files don't store port numbers, unfortunately, so
            # get rid of the port number if it's present.
            host = host.split(":")[0]

            mylog.debug("Looking for '%s %s' cookie in %s" % \
                  (host, path, self.cookie_file))
            self.cookie_jar.load(self.cookie_file, ignore_expires=True)

            try:
                cookie = self.cookie_jar._cookies[host][path]['rbsessionid']

                if not cookie.is_expired():
                    mylog.debug("Loaded valid cookie -- no login required")
                    return True

                mylog.debug("Cookie file loaded, but cookie has expired")
            except KeyError:
                mylog.debug("Cookie file loaded, but no cookie for this server")
        except IOError, error:
            mylog.debug("Couldn't load cookie file: %s" % error)

        return False

    def new_review_request(self, changenum, submit_as=None):
        """
        Creates a review request on a Review Board server, updating an
        existing one if the changeset number already exists.

        If submit_as is provided, the specified user name will be recorded as
        the submitter of the review request (given that the logged in user has
        the appropriate permissions).
        """
        try:
            mylog.debug("Attempting to create review request for %s" % changenum)
            data = { 'repository_path': self.info.path }

            if changenum:
                data['changenum'] = changenum

            if submit_as:
                mylog.debug("Submitting the review request as %s" % submit_as)
                data['submit_as'] = submit_as

            rsp = self.api_post('api/json/reviewrequests/new/', data)
        except APIError, e:
            rsp, = e.args

            if not self.options.diff_only:
                if rsp['err']['code'] == 204: # Change number in use
                    mylog.debug("Review request already exists. Updating it...")
                    rsp = self.api_post(
                        'api/json/reviewrequests/%s/update_from_changenum/' %
                        rsp['review_request']['id'])
                else:
                    raise e

        mylog.debug("Review request created")
        return rsp['review_request']

    def set_review_request_field(self, review_request, field, value):
        """
        Sets a field in a review request to the specified value.
        """
        rid = review_request['id']

        mylog.debug("Attempting to set field '%s' to '%s' for review request '%s'" %
              (field, value, rid))

        self.api_post('api/json/reviewrequests/%s/draft/set/' % rid, {
            field: value,
        })

    def get_review_request(self, rid):
        """
        Returns the review request with the specified ID.
        """
        rsp = self.api_get('api/json/reviewrequests/%s/' % rid)
        return rsp['review_request']

    def get_repositories(self):
        """
        Returns the list of repositories on this server.
        """
        rsp = self.api_get('/api/json/repositories/')
        return rsp['repositories']

    def get_repository_info(self, rid):
        """
        Returns detailed information about a specific repository.
        """
        rsp = self.api_get('/api/json/repositories/%s/info/' % rid)
        return rsp['info']

    def save_draft(self, review_request):
        """
        Saves a draft of a review request.
        """
        self.api_post("api/json/reviewrequests/%s/draft/save/" %
                      review_request['id'])
        mylog.debug("Review request draft saved")

    def upload_diff(self, review_request, diff_content, parent_diff_content):
        """
        Uploads a diff to a Review Board server.
        """
        mylog.debug("Uploading diff, size: %d" % len(diff_content))

        if parent_diff_content:
            mylog.debug("Uploading parent diff, size: %d" % len(parent_diff_content))

        fields = {}
        files = {}

        if self.info.base_path:
            fields['basedir'] = self.info.base_path

        files['path'] = {
            'filename': 'diff',
            'content': diff_content
        }

        if parent_diff_content:
            files['parent_diff_path'] = {
                'filename': 'parent_diff',
                'content': parent_diff_content
            }

        self.api_post('api/json/reviewrequests/%s/diff/new/' %
                      review_request['id'], fields, files)

    def publish(self, review_request):
        """
        Publishes a review request.
        """
        mylog.debug("Publishing")
        self.api_post('api/json/reviewrequests/%s/publish/' %
                      review_request['id'])

    def _get_server_info(self):
        if not self._server_info:
            self._server_info = self._info.find_server_repository_info(self)

        return self._server_info

    info = property(_get_server_info)

    def process_json(self, data):
        """
        Loads in a JSON file and returns the data if successful. On failure,
        APIError is raised.
        """
        rsp = json.loads(data)

        if rsp['stat'] == 'fail':
            raise APIError, rsp

        return rsp

    def http_get(self, path):
        """
        Performs an HTTP GET on the specified path, storing any cookies that
        were set.
        """
        mylog.debug('HTTP GETting %s' % path)

        url = self._make_url(path)

        try:
            rsp = urllib2.urlopen(url).read()
            self.cookie_jar.save(self.cookie_file)
            return rsp
        except urllib2.HTTPError, e:
            print "Unable to access %s (%s). The host path may be invalid" % \
                (url, e.code)
            try:
                mylog.debug(e.read())
            except AttributeError:
                pass
            die()

    def _make_url(self, path):
        """Given a path on the server returns a full http:// style url"""
        app = urlparse(self.url)[2]
        if path[0] == '/':
            url = urljoin(self.url, app[:-1] + path)
        else:
            url = urljoin(self.url, app + path)

        if not url.startswith('http'):
            url = 'http://%s' % url
        return url

    def api_get(self, path):
        """
        Performs an API call using HTTP GET at the specified path.
        """
        return self.process_json(self.http_get(path))

    def http_post(self, path, fields, files=None):
        """
        Performs an HTTP POST on the specified path, storing any cookies that
        were set.
        """
        if fields:
            debug_fields = fields.copy()
        else:
            debug_fields = {}

        if 'password' in debug_fields:
            debug_fields["password"] = "**************"
        url = self._make_url(path)
        mylog.debug('HTTP POSTing to %s: %s' % (url, debug_fields))

        content_type, body = self._encode_multipart_formdata(fields, files)
        headers = {
            'Content-Type': content_type,
            'Content-Length': str(len(body))
        }

        try:
            r = urllib2.Request(url, body, headers)
            data = urllib2.urlopen(r).read()
            self.cookie_jar.save(self.cookie_file)
            return data
        except urllib2.URLError, e:
            try:
                mylog.debug(e.read())
            except AttributeError:
                pass

            die("Unable to access %s. The host path may be invalid\n%s" % \
                (url, e))
        except urllib2.HTTPError, e:
            die("Unable to access %s (%s). The host path may be invalid\n%s" % \
                (url, e.code, e.read()))

    def api_post(self, path, fields=None, files=None):
        """
        Performs an API call using HTTP POST at the specified path.
        """
        return self.process_json(self.http_post(path, fields, files))

    def _encode_multipart_formdata(self, fields, files):
        """
        Encodes data for use in an HTTP POST.
        """
        BOUNDARY = mimetools.choose_boundary()
        content = ""

        fields = fields or {}
        files = files or {}

        for key in fields:
            content += "--" + BOUNDARY + "\r\n"
            content += "Content-Disposition: form-data; name=\"%s\"\r\n" % key
            content += "\r\n"
            content += fields[key] + "\r\n"

        for key in files:
            filename = files[key]['filename']
            value = files[key]['content']
            content += "--" + BOUNDARY + "\r\n"
            content += "Content-Disposition: form-data; name=\"%s\"; " % key
            content += "filename=\"%s\"\r\n" % filename
            content += "\r\n"
            content += value + "\r\n"

        content += "--" + BOUNDARY + "--\r\n"
        content += "\r\n"

        content_type = "multipart/form-data; boundary=%s" % BOUNDARY

        return content_type, content


class SCMClient(object):
    """
    A base representation of an SCM tool for fetching repository information
    and generating diffs.
    """
    def __init__(self, options):
        self.options = options

    def get_repository_info(self):
        return None

    def scan_for_server(self, repository_info):
        """
        Scans the current directory on up to find a .reviewboard file
        containing the server path.
        """
        raise RuntimeError
        server_url = self._get_server_from_config(user_config, repository_info)
        if server_url:
            return server_url

        for path in walk_parents(os.getcwd()):
            filename = os.path.join(path, ".reviewboardrc")
            if os.path.exists(filename):
                config = load_config_file(filename)
                server_url = self._get_server_from_config(config,
                                                          repository_info)
                if server_url:
                    return server_url

        return None

    def diff(self, args):
        """
        Returns the generated diff and optional parent diff for this
        repository.

        The returned tuple is (diff_string, parent_diff_string)
        """
        return (None, None)

    def diff_between_revisions(self, revision_range, args, repository_info):
        """
        Returns the generated diff between revisions in the repository.
        """
        return None

    def _get_server_from_config(self, config, repository_info):
        if 'REVIEWBOARD_URL' in config:
            return config['REVIEWBOARD_URL']
        elif 'TREES' in config:
            trees = config['TREES']
            if not isinstance(trees, dict):
                die("Warning: 'TREES' in config file is not a dict!")

            if repository_info.path in trees and \
               'REVIEWBOARD_URL' in trees[repository_info.path]:
                return trees[repository_info.path]['REVIEWBOARD_URL']

        return None

class SVNClient(SCMClient):
    """
    A wrapper around the svn Subversion tool that fetches repository
    information and generates compatible diffs.
    """
    def get_repository_info(self):
        if not check_install('svn help'):
            return None

        # Get the SVN repository path (either via a working copy or
        # a supplied URI)
        svn_info_params = ["svn", "info"]
        if self.options.repository_url:
            svn_info_params.append(self.options.repository_url)
        data = execute(svn_info_params,
                       ignore_errors=True)
        m = re.search(r'^Repository Root: (.+)$', data, re.M)
        if not m:
            return None

        path = m.group(1)

        m = re.search(r'^URL: (.+)$', data, re.M)
        if not m:
            return None

        base_path = m.group(1)[len(path):] or "/"

        m = re.search(r'^Repository UUID: (.+)$', data, re.M)
        if not m:
            return None

        return SvnRepositoryInfo(path, base_path, m.group(1))

    def scan_for_server(self, repository_info):
        # Scan first for dot files, since it's faster and will cover the
        # user's $HOME/.reviewboardrc
        server_url = super(SVNClient, self).scan_for_server(repository_info)
        if server_url:
            return server_url

        return self.scan_for_server_property(repository_info)

    def scan_for_server_property(self, repository_info):
        def get_url_prop(path):
            url = execute(["svn", "propget", "reviewboard:url", path]).strip()
            return url or None

        for path in walk_parents(os.getcwd()):
            if not os.path.exists(os.path.join(path, ".svn")):
                break

            prop = get_url_prop(path)
            if prop:
                return prop

        return get_url_prop(repository_info.path)

    def diff(self, files):
        """
        Performs a diff across all modified files in a Subversion repository.

        SVN repositories do not support branches of branches in a way that
        makes parent diffs possible, so we never return a parent diff
        (the second value in the tuple).
        """
        return (self.do_diff(["svn", "diff", "--diff-cmd=diff"] + files),
                None)

    def diff_between_revisions(self, revision_range, args, repository_info):
        """
        Performs a diff between 2 revisions of a Subversion repository.
        """
        if self.options.repository_url:
            revisions = revision_range.split(':')
            if len(revisions) < 1:
                return None
            elif len(revisions) == 1:
                revisions.append('HEAD')

            # if a new path was supplied at the command line, set it
            if len(args):
                repository_info.set_base_path(args[0])

            url = repository_info.path + repository_info.base_path

            old_url = url + '@' + revisions[0]
            new_url = url + '@' + revisions[1]

            return self.do_diff(["svn", "diff", "--diff-cmd=diff", old_url,
                                 new_url],
                                repository_info)
        # Otherwise, perform the revision range diff using a working copy
        else:
            return self.do_diff(["svn", "diff", "--diff-cmd=diff", "-r",
                                 revision_range],
                                repository_info)

    def do_diff(self, cmd, repository_info=None):
        """
        Performs the actual diff operation, handling renames and converting
        paths to absolute.
        """
        diff = execute(cmd, split_lines=True)
        diff = self.handle_renames(diff)
        diff = self.convert_to_absolute_paths(diff, repository_info)

        return ''.join(diff)

    def handle_renames(self, diff_content):
        """
        The output of svn diff is incorrect when the file in question came
        into being via svn mv/cp. Although the patch for these files are
        relative to its parent, the diff header doesn't reflect this.
        This function fixes the relevant section headers of the patch to
        portray this relationship.
        """

        # svn diff against a repository URL on two revisions appears to
        # handle moved files properly, so only adjust the diff file names
        # if they were created using a working copy.
        if self.options.repository_url:
            return diff_content

        result = []

        from_line = ""
        for line in diff_content:
            if line.startswith('--- '):
                from_line = line
                continue

            # This is where we decide how mangle the previous '--- '
            if line.startswith('+++ '):
                to_file, _ = self.parse_filename_header(line[4:])
                info       = self.svn_info(to_file)
                if info.has_key("Copied From URL"):
                    url       = info["Copied From URL"]
                    root      = info["Repository Root"]
                    from_file = urllib.unquote(url[len(root):])
                    result.append(from_line.replace(to_file, from_file))
                else:
                    result.append(from_line) #as is, no copy performed

            # We only mangle '---' lines. All others get added straight to
            # the output.
            result.append(line)

        return result


    def convert_to_absolute_paths(self, diff_content, repository_info):
        """
        Converts relative paths in a diff output to absolute paths.
        This handles paths that have been svn switched to other parts of the
        repository.
        """

        result = []

        for line in diff_content:
            front = None
            if line.startswith('+++ ') or line.startswith('--- ') or line.startswith('Index: '):
                front, line = line.split(" ", 1)

            if front:
                if line.startswith('/'): #already absolute
                    line = front + " " + line
                else:
                    # filename and rest of line (usually the revision
                    # component)
                    file, rest = self.parse_filename_header(line)

                    # If working with a diff generated outside of a working
                    # copy, then file paths are already absolute, so just
                    # add initial slash.
                    if self.options.repository_url:
                        path = urllib.unquote(
                            "%s/%s" % (repository_info.base_path, file))
                    else:
                        info = self.svn_info(file)
                        url  = info["URL"]
                        root = info["Repository Root"]
                        path = urllib.unquote(url[len(root):])

                    line = front + " " + path + rest

            result.append(line)

        return result

    def svn_info(self, path):
        """Return a dict which is the result of 'svn info' at a given path."""
        svninfo = {}
        for info in execute(["svn", "info", path],
                            split_lines=True):
            parts = info.strip().split(": ", 1)
            if len(parts) == 2:
                key, value = parts
                svninfo[key] = value

        return svninfo

    # Adapted from server code parser.py
    def parse_filename_header(self, s):
        parts = None
        if "\t" in s:
            # There's a \t separating the filename and info. This is the
            # best case scenario, since it allows for filenames with spaces
            # without much work.
            parts = s.split("\t")

        # There's spaces being used to separate the filename and info.
        # This is technically wrong, so all we can do is assume that
        # 1) the filename won't have multiple consecutive spaces, and
        # 2) there's at least 2 spaces separating the filename and info.
        if "  " in s:
            parts = re.split(r"  +", s)

        if parts:
            parts[1] = '\t' + parts[1]
            return parts

        # strip off ending newline, and return it as the second component
        return [s.split('\n')[0], '\n']

class MercurialClient(SCMClient):
    """
    A wrapper around the hg Mercurial tool that fetches repository
    information and generates compatible diffs.
    """
    def get_repository_info(self):
        if not check_install('hg --help'):
            return None

        data = execute(["hg", "root"], ignore_errors=True)
        if data.startswith('abort:'):
            # hg aborted => no mercurial repository here.
            return None

        # Elsewhere, hg root output give us the repository path.

        # We save data here to use it as a fallback. See below
        local_data = data.strip()

        svn = execute(["hg", "svn", "info", ], ignore_errors=True)

        if (not svn.startswith('abort:') and
            not svn.startswith("hg: unknown command")):
            self.type = 'svn'
            m = re.search(r'^Repository Root: (.+)$', svn, re.M)

            if not m:
                return None

            path = m.group(1)
            m2 = re.match(r'^(svn\+ssh|http|https)://([-a-zA-Z0-9.]*@)(.*)$',
                          path)
            if m2:
                path = '%s://%s' % (m2.group(1), m2.group(3))

            m = re.search(r'^URL: (.+)$', svn, re.M)

            if not m:
                return None

            base_path = m.group(1)[len(path):] or "/"
            return RepositoryInfo(path=path,
                                  base_path=base_path,
                                  supports_parent_diffs=True)

        self.type = 'hg'

        # We are going to search .hg/hgrc for the default path.
        file_name = os.path.join(local_data,'.hg', 'hgrc')

        if not os.path.exists(file_name):
            return RepositoryInfo(path=local_data, base_path='/',
                                  supports_parent_diffs=True)

        f = open(file_name)
        data = f.read()
        f.close()

        m = re.search(r'^default\s+=\s+(.+)$', data, re.M)

        if not m:
            # Return the local path, if no default value is found.
            return RepositoryInfo(path=local_data, base_path='/',
                                  supports_parent_diffs=True)

        path = m.group(1).strip()

        return RepositoryInfo(path=path, base_path='',
                              supports_parent_diffs=True)

    def diff(self, files):
        """
        Performs a diff across all modified files in a Mercurial repository.
        """
        # We don't support parent diffs with Mercurial yet, so return None
        # for the parent diff.
        if self.type == 'svn':
            return (execute(["hg", "svn", "diff", ]), None)

        return (execute(["hg", "diff"] + files), None)

    def diff_between_revisions(self, revision_range, args, repository_info):
        """
        Performs a diff between 2 revisions of a Mercurial repository.
        """
        if self.type != 'hg':
            raise NotImplementedError

        r1, r2 = revision_range.split(':')
        return execute(["hg", "diff", "-r", r1, "-r", r2])

def make_tempfile():
    """
    Creates a temporary file and returns the path. The path is stored
    in an array for later cleanup.
    """
    fd, tmpfile = mkstemp()
    os.close(fd)
    tempfiles.append(tmpfile)
    return tmpfile


def check_install(command):
    """
    Try executing an external command and return a boolean indicating whether
    that command is installed or not.  The 'command' argument should be
    something that executes quickly, without hitting the network (for
    instance, 'svn help' or 'git --version').
    """
    try:
        p = subprocess.Popen(command.split(' '),
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        return True
    except OSError:
        return False


def execute(command, env=None, split_lines=False, ignore_errors=False,
            extra_ignore_errors=(), translate_newlines=True):
    """
    Utility function to execute a command and return the output.
    """
    if isinstance(command, list):
        mylog.debug(subprocess.list2cmdline(command))
    else:
        mylog.debug(command)

    if env:
        env.update(os.environ)
    else:
        env = os.environ.copy()

    env['LC_ALL'] = 'en_US.UTF-8'
    env['LANGUAGE'] = 'en_US.UTF-8'

    if sys.platform.startswith('win'):
        p = subprocess.Popen(command,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=False,
                             universal_newlines=translate_newlines,
                             env=env)
    else:
        p = subprocess.Popen(command,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=False,
                             close_fds=True,
                             universal_newlines=translate_newlines,
                             env=env)
    if split_lines:
        data = p.stdout.readlines()
    else:
        data = p.stdout.read()
    rc = p.wait()
    if rc and not ignore_errors and rc not in extra_ignore_errors:
        die('Failed to execute command: %s\n%s' % (command, data))

    return data


def die(msg=None):
    """
    Cleanly exits the program with an error message. Erases all remaining
    temporary files.
    """
    for tmpfile in tempfiles:
        try:
            os.unlink(tmpfile)
        except:
            pass

    if msg:
        print msg

    sys.exit(1)


def walk_parents(path):
    """
    Walks up the tree to the root directory.
    """
    while os.path.splitdrive(path)[1] != os.sep:
        yield path
        path = os.path.dirname(path)


def load_config_file(filename):
    """
    Loads data from a config file.
    """
    config = {
        'TREES': {},
    }

    if os.path.exists(filename):
        try:
            execfile(filename, config)
        except:
            pass

    return config


def tempt_fate(server, tool, changenum, options, diff_content=None,
               parent_diff_content=None, submit_as=None, retries=3):
    """
    Attempts to create a review request on a Review Board server and upload
    a diff. On success, the review request path is displayed.
    """
    try:
        save_draft = False

        if options.rid:
            review_request = server.get_review_request(options.rid)
        else:
            review_request = server.new_review_request(changenum, submit_as)

        if options.target_groups:
            server.set_review_request_field(review_request, 'target_groups',
                                            options.target_groups)
            save_draft = True

        if options.target_people:
            server.set_review_request_field(review_request, 'target_people',
                                            options.target_people)
            save_draft = True

        if options.summary:
            server.set_review_request_field(review_request, 'summary',
                                            options.summary)
            save_draft = True

        if options.branch:
            server.set_review_request_field(review_request, 'branch',
                                            options.branch)
            save_draft = True

        if options.bugs_closed:
            server.set_review_request_field(review_request, 'bugs_closed',
                                            options.bugs_closed)
            save_draft = True

        if options.description:
            server.set_review_request_field(review_request, 'description',
                                            options.description)
            save_draft = True

        if options.testing_done:
            server.set_review_request_field(review_request, 'testing_done',
                                            options.testing_done)
            save_draft = True

        if save_draft:
            server.save_draft(review_request)
    except APIError, e:
        rsp, = e.args
        if rsp['err']['code'] == 103: # Not logged in
            retries = retries - 1

            # We had an odd issue where the server ended up a couple of
            # years in the future. Login succeeds but the cookie date was
            # "odd" so use of the cookie appeared to fail and eventually
            # ended up at max recursion depth :-(. Check for a maximum
            # number of retries.
            if retries >= 0:
                server.login(force=True)
                tempt_fate(server, tool, changenum, diff_content,
                           parent_diff_content, submit_as, retries=retries)
                return

        if options.rid:
            die("Error getting review request %s: %s (code %s)" % \
                (options.rid, rsp['err']['msg'], rsp['err']['code']))
        else:
            die("Error creating review request: %s (code %s)" % \
                (rsp['err']['msg'], rsp['err']['code']))


    if not server.info.supports_changesets or not options.change_only:
        try:
            server.upload_diff(review_request, diff_content,
                               parent_diff_content)
        except APIError, e:
            rsp, = e.args
            print "Error uploading diff: %s (%s)" % (rsp['err']['msg'],
                                                     rsp['err']['code'])
            mylog.debug(rsp)
            die("Your review request still exists, but the diff is not " +
                "attached.")

    if options.publish:
        server.publish(review_request)

    request_url = 'r/' + str(review_request['id'])
    review_url = urljoin(server.url, request_url)

    if not review_url.startswith('http'):
        review_url = 'http://%s' % review_url

    print "Review request #%s posted." % (review_request['id'],)
    print
    print review_url
    if not options.publish:
        print "This review will be inaccessible until you publish it!"

    return review_url


def parse_options(args):
    parser = OptionParser(usage="%prog [-pond] [-r review_id] [changenum]",
                          version="%prog " + VERSION)


    options, args = parser.parse_args(args)

    if options.description and options.description_file:
        sys.stderr.write("The --description and --description-file options "
                         "are mutually exclusive.\n")
        sys.exit(1)

    if options.description_file:
        if os.path.exists(options.description_file):
            fp = open(options.description_file, "r")
            options.description = fp.read()
            fp.close()
        else:
            sys.stderr.write("The description file %s does not exist.\n" %
                             options.description_file)
            sys.exit(1)

    if options.testing_done and options.testing_file:
        sys.stderr.write("The --testing-done and --testing-done-file options "
                         "are mutually exclusive.\n")
        sys.exit(1)

    if options.testing_file:
        if os.path.exists(options.testing_file):
            fp = open(options.testing_file, "r")
            options.testing_done = fp.read()
            fp.close()
        else:
            sys.stderr.write("The testing file %s does not exist.\n" %
                             options.testing_file)
            sys.exit(1)

    if options.repository_url and not options.revision_range:
        sys.stderr.write("The --repository-url option requires the "
                         "--revision-range option.\n")
        sys.exit(1)

    return args

def determine_client(options):

    repository_info = None
    tool = None

    # Try to find the SCM Client we're going to be working with.
    for tool in (SVNClient(options), MercurialClient(options)):
        repository_info = tool.get_repository_info()

        if repository_info:
            break

    if not repository_info:
        if options.repository_url:
            print "No supported repository could be access at the supplied url."
        else:
            print "The current directory does not contain a checkout from a"
            print "supported source code repository."
        sys.exit(1)

    # Verify that options specific to an SCM Client have not been mis-used.
    if options.change_only and not repository_info.supports_changesets:
        sys.stderr.write("The --change-only option is not valid for the "
                         "current SCM client.\n")
        sys.exit(1)

    if options.parent_branch and not repository_info.supports_parent_diffs:
        sys.stderr.write("The --parent option is not valid for the "
                         "current SCM client.\n")
        sys.exit(1)

    return (repository_info, tool)

def main(args, options):
    if 'USERPROFILE' in os.environ:
        homepath = os.path.join(os.environ["USERPROFILE"], "Local Settings",
                                "Application Data")
    elif 'HOME' in os.environ:
        homepath = os.environ["HOME"]
    else:
        homepath = ''

    # We don't want to do this...
    #globals()['user_config'] = \
    #    load_config_file(os.path.join(homepath, ".reviewboardrc"))

    # Load cookies
    cookie_file = os.path.join(homepath, ".post-review-cookies.txt")

    repository_info, tool = determine_client(options)

    server = ReviewBoardServer("http://review.enzotools.org/", repository_info,
                               cookie_file, options)

    if repository_info.supports_changesets:
        changenum = tool.get_changenum(args)
    else:
        changenum = None

    if options.revision_range:
        diff = tool.diff_between_revisions(options.revision_range, args,
                                           repository_info)
        parent_diff = None
    elif options.label and isinstance(tool, ClearCaseClient):
        diff, parent_diff = tool.diff_label(options.label)
    else:
        diff, parent_diff = tool.diff(args)

    if options.output_diff_only:
        print diff
        sys.exit(0)

    # Let's begin.
    server.login()

    review_url = tempt_fate(server, tool, changenum, options, diff_content=diff,
                            parent_diff_content=parent_diff,
                            submit_as=options.submit_as)


    # Load the review up in the browser if requested to:
    if options.open_browser:
        try:
            import webbrowser
            if 'open_new_tab' in dir(webbrowser):
                # open_new_tab is only in python 2.5+
                webbrowser.open_new_tab(review_url)
            elif 'open_new' in dir(webbrowser):
                webbrowser.open_new(review_url)
            else:
                os.system( 'start %s' % review_url )
        except:
            print 'Error opening review URL: %s' % review_url
