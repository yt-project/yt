from yt.funcs import issue_deprecation_warning

issue_deprecation_warning("The spectral_integrator module is deprecated. "
                          "'add_xray_emissivity_field' can now be imported "
                          "from the yt module.")

from yt.fields.xray_emission_fields import \
    add_xray_emissivity_field