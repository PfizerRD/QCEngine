%Mem={{ memory }}GB
%NProcShared={{ threads }}
%Chk=lig
# {{ theory }}/{{ basis }} {{ symmetry }} {{ driver }} {{ td_settings }} {%- for cmd in cmdline_extra %} {{ cmd }} {%- endfor %}

{{ title }}

{{ charge }} {{ multiplicity }}
{%- for element, coords in data %}
{{ element }}  {{ '{: .7f}'.format(coords[0]) }} {{ '{: .7f}'.format(coords[1]) }} {{ '{: .7f}'.format(coords[2]) }}
{%- endfor %}


{%- for cmd in add_input %}
{{ cmd }}
{%- endfor %}

