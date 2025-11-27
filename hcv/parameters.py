# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 17:33:20 2024

@author: kelly.maynard
"""

import numpy as np
from collections import OrderedDict
import pycountry

year_range = np.arange(2000, 2024)  # Years for data input
year_range_pop_data = np.arange(2000, 2024)

default_pops = {
    "PWID_males": {"label": "People Who Inject Drugs (males)", "type": "human"},
    "PWID_females": {"label": "People Who Inject Drugs (females)", "type": "human"},
    "0-9_males": {"label": "General Population 0-9 males", "type": "human"},
    "0-9_females": {"label": "General Population 0-9 females", "type": "human"},
    "10-17_males": {"label": "General Population 10-17 males", "type": "human"},
    "10-17_females": {"label": "General Population 10-17 females", "type": "human"},
    "18-64_males": {"label": "General Population 18-64 males", "type": "human"},
    "18-64_females": {"label": "General Population 18-64 females", "type": "human"},
    "65+_males": {"label": "General Population 65+ males", "type": "human"},
    "65+_females": {"label": "General Population 65+ females", "type": "human"},
    "Prisoners_males": {"label": "Prisoners (males)", "type": "human"},
    "Prisoners_females": {"label": "Prisoners (females)", "type": "human"},
}

general_populations = ["PWID", "0-9", "10-17", "18-64", "65+", "Prisoners"]

# Population calculations setup
sub_pops = [
    "Prisoners",
    "PWID",
]  # ,'Hospital','Indigenous'] #pops to be subtracted from gen pop size
default_pops_inter = [
    "Prisoners_males",
    "Prisoners_females",
]  # default populations with idu self interactions

# default transfers - really important for these all to be in the correct order (ie. age, idu, inc) because D.transfers is indexed by number rather than parameter
default_transfer_codes = OrderedDict()
default_transfer_codes = {
    "age": "Aging",
    "idu_status": "Injecting Drug Use Relapse or Cessation",
    "inc": "Incarceration",
}
default_transfers = OrderedDict()
default_transfers = {
    "age": [
        ("18-64_males", "65+_males"),
        ("18-64_females", "65+_females"),
        ("10-17_males", "18-64_males"),
        ("10-17_females", "18-64_females"),
        ("0-9_males", "10-17_males"),
        ("0-9_females", "10-17_females"),
    ],
    "idu_status": [
        ("PWID_males", "18-64_males"),
        ("18-64_males", "PWID_males"),
        ("PWID_females", "18-64_females"),
        ("18-64_females", "PWID_females"),
    ],
    "inc": [
        ("PWID_males", "Prisoners_males"),
        ("PWID_females", "Prisoners_females"),
        ("Prisoners_males", "PWID_males"),
        ("Prisoners_females", "PWID_females"),
        ("18-64_males", "Prisoners_males"),
        ("18-64_females", "Prisoners_females"),
        ("Prisoners_males", "18-64_males"),
        ("Prisoners_females", "18-64_females"),
    ],
}


# Country to ISO dictionary
replace_country_names = {
    "Antigua & Barbuda": "Antigua and Barbuda",
    "Sao Tome & Principe": "Sao Tome and Principe",
    "Republic of Korea": "Korea, Republic of",
    "Holy See": "Holy See (Vatican City State)",
    "Taiwan": "Taiwan, Province of China",
    "Bosnia & Herzegovina": "Bosnia and Herzegovina",
    "Bolivia": "Bolivia, Plurinational State of",
    "Venezuela": "Venezuela, Bolivarian Republic of",
    "Micronesia": "Micronesia, Federated States of",
    "Cote d'Ivoire": "Côte d'Ivoire",
    "Syria": "Syrian Arab Republic",
    "Czech Republic": "Czechia",
    "North Korea": "Korea, Democratic People's Republic of",
    "Laos": "Lao People's Democratic Republic",
    "United States of America": "United States",
    "State of Palestine": "Palestine, State of",
    "United Republic of Tanzania": "Tanzania, United Republic of",
    "St Vincent & Grenadines": "Saint Vincent and the Grenadines",
    "Iran": "Iran, Islamic Republic of",
    "Saint Kitts & Nevis": "Saint Kitts and Nevis",
    "Moldova": "Moldova, Republic of",
    "DR Congo": "Congo, The Democratic Republic of the",
    "Turkey": "Türkiye",
    "Trinidad & Tobago": "Trinidad and Tobago",
}  # NOTE: These keys are reversed

replace_country_names = {v: k for k, v in replace_country_names.items()}

# Generate a dictionary for country ISO codes
country_to_iso = {}

for country in pycountry.countries:
    country_to_iso[country.name] = country.alpha_3

country_to_iso["Kosovo"] = "XKX"

for country in replace_country_names.keys():
    country_to_iso[replace_country_names[country]] = country_to_iso.pop(country)

iso_to_country = {v: k for k, v in country_to_iso.items()}

country_to_iso["Russia"] = "RUS"
country_to_iso["United States"] = "USA"
country_to_iso["Vietnam"] = "VNM"
