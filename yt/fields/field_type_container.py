"""
A proxy object for field descriptors, usually living as ds.fields.
"""

import inspect
import sys
import textwrap
import weakref

from yt.fields.derived_field import DerivedField

if sys.version_info >= (3, 8):
    from functools import cached_property
else:
    from yt._maintenance.backports import cached_property


def _fill_values(values):
    value = (
        '<div class="rendered_html jp-RenderedHTMLCommon">'
        + "<table><thead><tr><th>Name</th><th>Type</th>"
        + "<th>Value</th></tr></thead><tr><td>"
        + "</td></tr><tr><td>".join(
            "{}</td><td>{}</td><td>{}".format(
                v, type(values[v]).__name__, str(values[v])
            )
            for v in sorted(values)
        )
        + "</td></tr></table></div>"
    )
    return value


class FieldTypeContainer:
    def __init__(self, ds):
        self.ds = weakref.proxy(ds)

    def __getattr__(self, attr):
        ds = self.__getattribute__("ds")
        fnc = FieldNameContainer(ds, attr)
        if len(dir(fnc)) == 0:
            return self.__getattribute__(attr)
        return fnc

    @cached_property
    def field_types(self):
        return {t for t, n in self.ds.field_info}

    def __dir__(self):
        return list(self.field_types)

    def __iter__(self):
        for ft in self.field_types:
            fnc = FieldNameContainer(self.ds, ft)
            if len(dir(fnc)) == 0:
                yield self.__getattribute__(ft)
            else:
                yield fnc

    def __contains__(self, obj):
        ob = None
        if isinstance(obj, FieldNameContainer):
            ob = obj.field_type
        elif isinstance(obj, str):
            ob = obj

        return ob in self.field_types

    def _ipython_display_(self):
        import ipywidgets
        from IPython.display import display

        fnames = []
        children = []
        for ftype in sorted(self.field_types):
            fnc = getattr(self, ftype)
            children.append(ipywidgets.Output())
            with children[-1]:
                display(fnc)
            fnames.append(ftype)
        tabs = ipywidgets.Tab(children=children)
        for i, n in enumerate(fnames):
            tabs.set_title(i, n)
        display(tabs)


class FieldNameContainer:
    def __init__(self, ds, field_type):
        self.ds = ds
        self.field_type = field_type

    def __getattr__(self, attr):
        ft = self.__getattribute__("field_type")
        ds = self.__getattribute__("ds")
        if (ft, attr) not in ds.field_info:
            return self.__getattribute__(attr)
        return ds.field_info[ft, attr]

    def __dir__(self):
        return [n for t, n in self.ds.field_info if t == self.field_type]

    def __iter__(self):
        for t, n in self.ds.field_info:
            if t == self.field_type:
                yield self.ds.field_info[t, n]

    def __contains__(self, obj):
        if isinstance(obj, DerivedField):
            if self.field_type == obj.name[0] and obj.name in self.ds.field_info:
                # e.g. from a completely different dataset
                if self.ds.field_info[obj.name] is not obj:
                    return False
                return True
        elif isinstance(obj, tuple):
            if self.field_type == obj[0] and obj in self.ds.field_info:
                return True
        elif isinstance(obj, str):
            if (self.field_type, obj) in self.ds.field_info:
                return True
        return False

    def _ipython_display_(self):
        import ipywidgets
        from IPython.display import Markdown, display

        names = dir(self)
        names.sort()

        def change_field(_ftype, _box, _var_window):
            def _change_field(event):
                fobj = getattr(_ftype, event["new"])
                _box.clear_output()
                with _box:
                    display(
                        Markdown(
                            data="```python\n"
                            + textwrap.dedent(fobj.get_source())
                            + "\n```"
                        )
                    )
                values = inspect.getclosurevars(fobj._function).nonlocals
                _var_window.value = _fill_values(values)

            return _change_field

        flist = ipywidgets.Select(options=names, layout=ipywidgets.Layout(height="95%"))
        source = ipywidgets.Output(layout=ipywidgets.Layout(width="100%", height="9em"))
        var_window = ipywidgets.HTML(value="Empty")
        var_box = ipywidgets.Box(
            layout=ipywidgets.Layout(width="100%", height="100%", overflow_y="scroll")
        )
        var_box.children = [var_window]
        ftype_tabs = ipywidgets.Tab(
            children=[source, var_box],
            layout=ipywidgets.Layout(flex="2 1 auto", width="auto", height="95%"),
        )
        ftype_tabs.set_title(0, "Source")
        ftype_tabs.set_title(1, "Variables")
        flist.observe(change_field(self, source, var_window), "value")
        display(
            ipywidgets.HBox(
                [flist, ftype_tabs], layout=ipywidgets.Layout(height="14em")
            )
        )
