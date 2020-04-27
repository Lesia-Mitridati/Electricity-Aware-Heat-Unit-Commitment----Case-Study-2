# Electricity-Aware-Heat-Unit-Commitment----Online-Appendix
Online appendix to companion paper submitted to IJOO, Danish Energy system case study setup, and Python code for EAHUC model. 

Nomenclature:

hiy: hour in year (without daylight saving time change)
diy: day in year
hid: hour in day
Pmin/Pmax: minimum and maximum electricity output (MW)
Qmin/Qmax: minimum and maximum heat output (MW)
Ce: electricity (linear) cost parameter (EUR/MWh)
Ch: heat (linear) cost parameter (EUR/MWh)
C0: No-load costs (EUR)
Cs: start-up costs (EUR)
Ns: number of hours with distinct start-up costs (h)
Emax: maximum energy stored (MWh)
Relec: electricity efficiency of fuel of CHPs (-)
Rheat: heat efficiency of fuel of CHPs (-)
Fmin/Fmax: minimum and maximum fuel intake (MW)
Rmin: minimum power/heat ratio of CHPs or inverse coefficient of performance of heat pumps, s.t. P = Rmin*Q (-)

Electricity market zones:

DK1-DK2: Denmark
SE1-SE4: Sweden
NO1-NO2: Norway
SE: Germany

Heat nodes:

TVIS_1-TVIS_6: TVIS district heating network
KBH1-KBH_14: Greater Copenhagen district heating network
AAR_1-AAR_9: Aarhus district heating network

Input data files:

ATC_upper.csv: Upper bound on ATC between market zones for each hiy
ATC_lower.csv: Lower bound on ATC between market zones for each hiy
elec_load.csv: Electricity consumption in Denmark (market zones DK1 and DK2) for each hiy
heat_load.csv: Heat consumption in nodes for each hiy
elec_network_topology.csv: 0/1 table of connections between market zones (1=connected, 0=disconnected)
heat_network_topology.csv: 0/1 table of connections between nodes (1=connected, 0=disconnected)
solar_prod.csv: Electricity production from PV in Denmark (market zones DK1 and DK2) for each hiy
wind_prod.csv: Electricity production from wind in Denmark (market zones DK1 and DK2) for each hiy
NE_data_heat_elec.xslx: Excel file of all input data

Python code: 

build_bids.py: estimate electricity prices and compute heat bids (price-quantity) of heat pumps and CHPs, for each hour of the following day
HUC_no_time_delay.py: electricity-aware, combined, and decoupled heat unit commitment and heat market clearings, with time delays in pipelines fixed to 0, for each hour of the following day
HUC.py: electricity-aware, combined and decoupled heat unit commitment and heat market clearings, with variable time delays in pipelines, for each hour of the following day
EM.py: electricity market clearing for eachhour of the following day

