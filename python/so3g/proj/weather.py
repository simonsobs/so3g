import ephem


class Weather(dict):
    """This is a thin wrapper around a dict or some other convenient
    container, to simplify passing telescope atmospheric weather
    information to various consumers.

    Important keys:

    - 'temperature': ambient temperature, in degrees Celsius.
    - 'pressure': air pressure, in mBar.
    - 'humidity': relative humidity, as fraction of 1.

    If you find that some consumer prefers some other units (K and
    atm, for example...), then define a function that makes the
    necessary conversions.

    """

    default = {'temperature': 0.,
               'pressure': 0.,
               'humidity': 0.}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for k, v in self.default.items():
            if k not in self:
                self[k] = v

    def apply(self, something):
        """Helper method to load configure an object with this weather info.
        Behavior depends on the type of observer.  The types below
        will be checked, in order.

        - ephem.Observer -- sets pressure and temperature.

        """
        if isinstance(something, ephem.Observer):
            something.pressure = self['pressure']
            something.temp = self['temperature']
        else:
            raise TypeError('No Weather application procedure defined for '
                            'objects of type %s' % something.__class__)

    def to_qpoint(self):
        """Extract a dict appropriate for passing to qpoint functions.

        """
        return {k: self[k] for k in ['temperature', 'pressure', 'humidity']}


def weather_factory(settings):
    """
    Helper function that returns some useful Weather objects.
    """
    if settings == 'vacuum':
        return Weather()
    elif settings in ['toco', 'act', 'so', 'sa']:
        return Weather({'temperature': 0.,
                        'humidity': 0.2,
                        'pressure': 550.})
    else:
        raise ValueError('No factory machinery for settings=%s' %
                         str(settings))
