import os
import hglib
import pkg_resources
from email.utils import parseaddr

templates = {"header": r"""<Project xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#" xmlns="http://usefulinc.com/ns/doap#" xmlns:foaf="http://xmlns.com/foaf/0.1/" xmlns:admin="http://webns.net/mvcb/">
 <name>The yt Project</name>
 <shortname>yt</shortname>
 <shortdesc>Multi-resolution volumetric analysis</shortdesc>
 <description>yt is a python package for analyzing and visualizing volumetric, multi-resolution data.  Originally developed for astrophysical simulations, it is flexible enough to work with data from other domains such as weather, nuclear engineering, seismology and molecular dynamics.</description>
 <homepage rdf:resource="http://yt-project.org/" />
 <download-page rdf:resource="https://pypi.python.org/pypi/yt/" />
 <download-mirror rdf:resource="http://bitbucket.org/yt_analysis/yt/" />
 <bug-database rdf:resource="http://bitbucket.org/yt_analysis/yt/issues" />
 <programming-language>python</programming-language>
 <programming-language>cython</programming-language>
 <license rdf:resource="http://usefulinc.com/doap/licenses/bsd" />
 """,
 "foaf": r"""<foaf:Person>
     <foaf:name>%(realname)s</foaf:name>
  </foaf:Person>
  """,
  "release": r"""
	<release>
		<Version>
			<name>%(name)s</name>
			<created>%(date)s</created>
			<revision>%(revision)s</revision>
		</Version>
	</release>
   """,
  "footer": r"""
 <repository> 
   <HgRepository>
     <browse rdf:resource='https://bitbucket.org/yt_analysis/yt/src' />
     <location rdf:resource='https://bitbucket.org/yt_analysis/yt' />
   </HgRepository>
 </repository> 
</Project>
"""
}

releases = [
    ("1.0.1", "2008-10-25"),
    ("1.5"  , "2009-11-04"),
    ("1.6"  , "2010-01-22"),
    ("1.6.1", "2010-02-11"),
    ("1.7"  , "2010-06-27"),
    ("2.0"  , "2011-01-17"),
    ("2.0.1", "2011-01-20"),
    ("2.1"  , "2011-04-06"),
    ("2.2"  , "2011-09-02"),
    ("2.3"  , "2011-12-15"),
    ("2.4"  , "2012-08-02"),
    ("2.5"  , "2013-03-01"),
    ("2.5.1", "2013-03-31"),
    ("2.5.2", "2013-05-01"),
    ("2.5.3", "2013-06-03"),
    ("2.5.4", "2013-07-02"),
    ("2.5.5", "2013-08-23"),
    ("2.6"  , "2013-11-23"),
    ("2.6.1", "2013-12-03"),
    ("2.6.2", "2014-02-28"),
    ("2.6.3", "2014-07-23"),
    ("3.0"  , "2014-08-04"),
    ("3.0.1", "2014-09-01"),
    ("3.0.2", "2014-10-03"),
    ("3.0.3", "2014-11-03"),
    ("3.1"  , "2015-01-14"),
    ("3.2"  , "2015-07-24"),
    ("3.2.1", "2015-09-09"),
    ("3.2.2", "2015-11-13"),
]

yt_provider = pkg_resources.get_provider("yt")
yt_path = os.path.dirname(yt_provider.module_path)

name_mappings = {
        # Sometimes things get filtered out by hgchurn pointing elsewhere.
        # So we can add them back in manually
        "andrew.wetzel@yale.edu" : "Andrew Wetzel",
        "df11c@my.fsu.edu": "Daniel Fenn",
        "dnarayan@haverford.edu": "Desika Narayanan",
        "jmtomlinson95@gmail.com": "Joseph Tomlinson",
        "kaylea.nelson@yale.edu": "Kaylea Nelson",
        "tabel@slac.stanford.edu": "Tom Abel",
        "pshriwise": "Patrick Shriwise",
        "jnaiman": "Jill Naiman",
        "gsiisg": "Geoffrey So",
        "dcollins4096@gmail.com": "David Collins",
        "bcrosby": "Brian Crosby",
        "astrugarek": "Antoine Strugarek",
        "AJ": "Allyson Julian",
}

name_ignores = ["convert-repo"]

lastname_sort = lambda a: a.rsplit(None, 1)[-1]

def developer_names():
    cmd = hglib.util.cmdbuilder("churn", "-c")
    c = hglib.open(yt_path)
    emails = set([])
    for dev in c.rawcommand(cmd).split("\n"):
        if len(dev.strip()) == 0: continue
        emails.add(dev.rsplit(None, 2)[0])
    print("Generating real names for {0} emails".format(len(emails)))
    names = set([])
    for email in sorted(emails):
        if email in name_ignores:
            continue
        if email in name_mappings:
            names.add(name_mappings[email])
            continue
        cset = c.log(revrange="last(author('%s'))" % email)
        if len(cset) == 0:
            print("Error finding {0}".format(email))
            realname = email
        else:
            realname, addr = parseaddr(cset[0][4])
        if realname == '':
            realname = email
        if realname in name_mappings:
            names.add(name_mappings[realname])
            continue
        realname = realname.decode('utf-8')
        realname = realname.encode('ascii', 'xmlcharrefreplace')
        names.add(realname)
    #with open("devs.txt", "w") as f:
    #    for name in sorted(names, key=lastname_sort):
    #        f.write("%s\n" % name)
    devs = list(names)
    devs.sort(key=lastname_sort)
    return devs

def generate_doap():
    dev_names = developer_names()
    with open("doap.xml", "w") as f:
        f.write(templates["header"])
        for dev_name in dev_names:
            f.write("<developer>\n")
            f.write(templates["foaf"] % {'realname': dev_name})
            f.write("</developer>\n")
        for release in releases:
            f.write(templates["release"] % {
                'name': "yt " + release[0], 'revision': release[0], 'date': release[1]}
            )
        f.write(templates["footer"])

if __name__ == "__main__":
    generate_doap()
