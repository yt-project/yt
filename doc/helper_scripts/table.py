contents = [
    (
        "Getting Started",
        [
            ("welcome/index.html", "Welcome to yt!", "What's yt all about?"),
            (
                "orientation/index.html",
                "yt Orientation",
                "Quickly get up and running with yt: zero to sixty.",
            ),
            (
                "help/index.html",
                "How to Ask for Help",
                "Some guidelines on how and where to ask for help with yt",
            ),
            (
                "workshop.html",
                "Workshop Tutorials",
                "Videos, slides and scripts from the 2012 workshop covering many "
                + "aspects of yt, from beginning to advanced.",
            ),
        ],
    ),
    (
        "Everyday yt",
        [
            (
                "analyzing/index.html",
                "Analyzing Data",
                "An overview of different ways to handle and process data: loading "
                + "data, using and manipulating objects and fields, examining and "
                + "manipulating particles, derived fields, generating processed data, "
                + "time series analysis.",
            ),
            (
                "visualizing/index.html",
                "Visualizing Data",
                "An overview of different ways to visualize data: making projections, "
                + "slices, phase plots, streamlines, and volume rendering; modifying "
                + "plots; the fixed resolution buffer.",
            ),
            (
                "interacting/index.html",
                "Interacting with yt",
                "Different ways -- scripting, GUIs, prompts, explorers -- to explore "
                + "your data.",
            ),
        ],
    ),
    (
        "Advanced Usage",
        [
            (
                "advanced/index.html",
                "Advanced yt usage",
                "Advanced topics: parallelism, debugging, ways to drive yt, "
                + "developing",
            ),
            (
                "getting_involved/index.html",
                "Getting Involved",
                "How to participate in the community, contribute code and share "
                + "scripts",
            ),
        ],
    ),
    (
        "Reference Materials",
        [
            (
                "cookbook/index.html",
                "The Cookbook",
                "A bunch of illustrated examples of how to do things",
            ),
            (
                "reference/index.html",
                "Reference Materials",
                "A list of all bundled fields, API documentation, the Change Log...",
            ),
            ("faq/index.html", "FAQ", "Frequently Asked Questions: answered for you!"),
        ],
    ),
]

heading_template = r"""
<h2>%s</h2>
<table class="contentstable" align="center">
%s
</table>
"""

subheading_template = r"""
  <tr valign="top">
    <td width="25%%">
      <p>
        <a href="%s">%s</a>
      </p>
    </td>
    <td width="75%%">
      <p class="linkdescr">%s</p>
    </td>
  </tr>
"""

t = ""
for heading, items in contents:
    s = ""
    for subheading in items:
        s += subheading_template % subheading
    t += heading_template % (heading, s)
print(t)
