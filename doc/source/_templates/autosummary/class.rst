{% extends "!autosummary/class.rst" %}

{% block methods %}
{% if methods %}
.. autosummary::
   :toctree:
   {% for item in methods %}
       ~{{ name }}.{{ item }}
   {%- endfor %}

{% endif %}
{% endblock %}

{% block attributes %}
{% if attributes %}
.. autosummary::
   :toctree:
   {% for item in attributes %}
       ~{{ name }}.{{ item }}
   {%- endfor %}

{% endif %}
{% endblock %}
