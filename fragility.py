#!/usr/bin/env python3

'''
Classes for handling of the fragility functions
and the damage states.
'''

import collections
import json
import re

from scipy.stats import lognorm
import numpy as np


class LogncdfFactory():
    '''
    This is function factory for the log normal cdf.
    '''

    def __call__(self, mean, stddev):
        func = lognorm(scale=np.exp(mean), s=stddev)
        return func.cdf


SUPPORTED_FRAGILITY_FUNCTION_FACTORIES = {
    'logncdf': LogncdfFactory(),
}


class DamageState():
    '''
    Class to represent the damage states.
    '''
    def __init__(self,
                 taxonomy,
                 from_state,
                 to_state,
                 intensity_field,
                 intensity_unit,
                 fragility_function):
        self.taxonomy = taxonomy
        self.from_state = from_state
        self.to_state = to_state
        self.intensity_field = intensity_field
        self.intensity_unit = intensity_unit

        self.fragility_function = fragility_function

    def get_probability_for_intensity(self, intensity, units):
        '''
        Returns the probabilit value for the given
        intensity.

        The intensity and units are given as dicts, for example:
        intensity = {
            'PGA': 1.0,
            'STDDEV_PGA': 7.0
        }
        units = {
            'PGA': 'g',
            'STDDEV_PGA': 'g'
        }

        This method throws an exception if the unit for the
        fragility function is not the expected one.
        '''
        field = self.intensity_field.upper()
        value = intensity[field]
        unit = units[field]

        if unit != self.intensity_unit:
            raise Exception('Not supported unit')

        return self.fragility_function(value)


class Fragility():
    '''
    Class to represent all of the fragility data.
    '''

    def __init__(self, data):
        self._data = data

    @classmethod
    def from_file(cls, json_file):
        '''
        Reads the data from a given json file.
        '''
        with open(json_file, 'rt') as input_file:
            data = json.load(input_file)
        return cls(data)

    def to_fragility_provider(self):
        '''
        Transforms the data, so that a
        provider for the supported taxonomies
        and the damage states (with the fragility functions)
        are returned.
        '''
        damage_states_by_taxonomy = collections.defaultdict(list)

        shape = self._data['meta']['shape']

        for dataset in self._data['data']:
            taxonomy = dataset['taxonomy']
            intensity_field = dataset['imt']
            intensity_unit = dataset['imu']
            for damage_state_mean_key in [
                    k for k in dataset.keys()
                    if k.startswith('D')
                    and k.endswith('_mean')]:
                #
                # the data is in the format
                # D1_mean, D2_mean, D3_mean
                # (as there are is no from data state at the moment)
                # but this code can also handle them the in the way
                # D01, so that it is the damage state from 0 to 1 or
                # D_0_1 or D0_1
                #
                to_state = int(
                    re.search(
                        r'(\d)_mean$',
                        damage_state_mean_key
                    ).group(1))
                from_state = int(
                    re.search(
                        r'^D_?(\d)_',
                        damage_state_mean_key
                    ).group(1))

                if to_state == from_state:
                    # there is no from state given
                    # both regexp read the same value
                    from_state = 0

                mean = dataset[damage_state_mean_key]
                stddev_key = damage_state_mean_key.replace('_mean', '_stddev')
                stddev = dataset[stddev_key]

                damage_state = DamageState(
                    taxonomy=taxonomy,
                    from_state=from_state,
                    to_state=to_state,
                    intensity_field=intensity_field,
                    intensity_unit=intensity_unit,
                    fragility_function=SUPPORTED_FRAGILITY_FUNCTION_FACTORIES[
                        shape](mean, stddev)
                )

                damage_states_by_taxonomy[taxonomy].append(damage_state)

        schema = self._data['meta']['id']

        return FragilityProvider(damage_states_by_taxonomy, schema)

class FragilityProvider():
    '''
    Class to give access to the taxonomies and
    the damage states with the fragility functions.
    '''
    def __init__(self, damage_states_by_taxonomy, schema):
        self._damage_states_by_taxonomy = damage_states_by_taxonomy
        self._schema = schema

    def get_damage_states_for_taxonomy(self, taxonomy):
        '''
        Returns all the damage states for the given
        taxonomy.
        '''
        return self._damage_states_by_taxonomy[taxonomy]

    def get_taxonomies(self):
        '''
        Returns the taxonomies from the data.
        '''
        return self._damage_states_by_taxonomy.keys()

    def get_schema(self):
        return self._schema
