import yt
from yt.analysis_modules.cosmological_observation.light_ray.api import LightRay
from yt.analysis_modules.absorption_spectrum.api import AbsorptionSpectrum
from yt.analysis_modules.absorption_spectrum.api import generate_total_fit

# Define a field to simulate OVI based on a constant relationship to HI
# Do *NOT* use this for science, because this is not how OVI actually behaves;
# it is just an example.

def _OVI_number_density(field, data):
    return data['H_number_density']*2.0

# Define a function that will accept a ds and add the new field
# defined above.  This will be given to the LightRay below.
def setup_ds(ds):
    ds.add_field(("gas","O_p5_number_density"),
                 function=_OVI_number_density,
                 units="cm**-3")

# Define species and associated parameters to add to continuum
# Parameters used for both adding the transition to the spectrum
# and for fitting
# Note that for single species that produce multiple lines
# (as in the OVI doublet), 'numLines' will be equal to the number
# of lines, and f,gamma, and wavelength will have multiple values.

HI_parameters = {'name': 'HI',
                 'field': 'H_number_density',
                 'f': [.4164],
                 'Gamma': [6.265E8],
                 'wavelength': [1215.67],
                 'mass': 1.00794,
                 'numLines': 1,
                 'maxN': 1E22, 'minN': 1E11,
                 'maxb': 300, 'minb': 1,
                 'maxz': 6, 'minz': 0,
                 'init_b': 30,
                 'init_N': 1E14}

OVI_parameters = {'name': 'OVI',
                  'field': 'O_p5_number_density',
                  'f': [.1325, .06580],
                  'Gamma': [4.148E8, 4.076E8],
                  'wavelength': [1031.9261, 1037.6167],
                  'mass': 15.9994,
                  'numLines': 2,
                  'maxN': 1E17, 'minN': 1E11,
                  'maxb': 300, 'minb': 1,
                  'maxz': 6, 'minz': 0,
                  'init_b': 20,
                  'init_N': 1E12}

species_dicts = {'HI': HI_parameters, 'OVI': OVI_parameters}

# Create a LightRay object extending from z = 0 to z = 0.1
# and use only the redshift dumps.
lr = LightRay('enzo_cosmology_plus/AMRCosmology.enzo',
              'Enzo', 0.0, 0.1,
              use_minimum_datasets=True,
              time_data=False
              )

# Get all fields that need to be added to the light ray
fields = ['temperature']
for s, params in species_dicts.items():
    fields.append(params['field'])

# Make a light ray, and set njobs to -1 to use one core
# per dataset.
lr.make_light_ray(seed=123456780,
                  solution_filename='lightraysolution.txt',
                  data_filename='lightray.h5',
                  fields=fields, setup_function=setup_ds,
                  njobs=-1)

# Create an AbsorptionSpectrum object extending from
# lambda = 900 to lambda = 1800, with 10000 pixels
sp = AbsorptionSpectrum(900.0, 1400.0, 50000)

# Iterate over species
for s, params in species_dicts.items():
    # Iterate over transitions for a single species
    for i in range(params['numLines']):
        # Add the lines to the spectrum
        sp.add_line(s, params['field'],
                    params['wavelength'][i], params['f'][i],
                    params['Gamma'][i], params['mass'],
                    label_threshold=1.e10)


# Make and save spectrum
wavelength, flux = sp.make_spectrum('lightray.h5',
                                    output_file='spectrum.h5',
                                    line_list_file='lines.txt',
                                    use_peculiar_velocity=True)


# Define order to fit species in
order_fits = ['OVI', 'HI']

# Fit spectrum and save fit
fitted_lines, fitted_flux = generate_total_fit(wavelength,
                                               flux, order_fits, species_dicts,
                                               output_file='spectrum_fit.h5')
