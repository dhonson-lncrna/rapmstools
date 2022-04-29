import pandas as pd

# List of RNA targets
targets = ('XIST', 
           '7SL',
           'RMRP',
           'U7',
           '7SK',
           'U2',
           'U6',
           'U1',
           'GO Splice')

# Lists of known interactors
# Xist, from McHugh et al, Nature 2015 (https://doi.org/10.1038/nature14443)
xist = ('Spen',
       'Rbm15',
       'Myef2',
       'Celf1',
       'Hnrnpc',
       'Lbr',
       'Hnrnpu',
       'Raly',
       'Hnrnpm',
       'Ptbp1')

# 7SL, from Luirink and Sinning, BBA-MCR 2004 (https://doi.org/10.1016/j.bbamcr.2004.03.013)
sevensl = ('Srp9',
           'Srp14',
           'Srp19'
           'Srp68',
           'Srp72',
           'Srp54',
           'Srpr')

# RMRP, from Jarrous, Trends in Genetics 2017 (https://doi.org/10.1016/j.tig.2017.06.006)
# I could not find a list of all the mouse genes, so these are human
rmrp = ('Pop1',
        'Pop4'
       'Rpp38',
       'Rpp29',
       'Rpp21',
       'Rpp40',
       'Rpp14',
       'Pop5',
       'Rpp30',
       'Rpp25',
       'Rpp20')

# U7, from Schumperli and Pillai, CMLS 2004 (https://doi.org/10.1007/s00018-004-4190-0)
# Lsm10 and Lsm11 are specific to U7
useven = ('Lsm10',
         'Lsm11',
         'Snrpf',
         'Snrpb',
         'Snrpd3',
         'Snrpg',
         'Snrpe',
         'Zfp100',
         'Slbp')

# 7SK, from Brogie and Price, Nucleic Acids Research 2017 (doi: 10.1093/nar/gkx262)
# Also from Barrandon et al, MCM 2007 (https://doi.org/10.1128/MCB.00975-07)
sevensk = ('Hexim1',
          'Mepce',
          'Larp7',
          'Hnrnpq',
          'Hnrnpr',
          'Hnrnpa1',
          'Hnrnpa2',
          'Cdk9',
          'Cyct1')

# U2, from Zhang et al, Nature 2020 (https://doi.org/10.1038/s41586-020-2344-3)
# Also Scofield and Lynch, Mol Bio Evol 2008 (doi:10.1093/molbev/msn175)
utwo = ('Snrpa1',
       'Snrpb2',
       'Sf3a3',
       'Sf3b3',
       'Sf3b2',
       'Sf3b1',
       'Prp5',
       'Sf3b6',
       'Sf3b5',
       'Tat-sf1',
       'Sf3b4',
       'Sf3a2',
       'Sf3a1',
       'Snrpb',
       'Snrpd1',
       'Snrpd2',
       'Snrpd3'
       'Snrpg',
       'Snrpe',
       'Snrpf',
       'Snrpg')

# U6, Montemayor et al, Nature Comms 2018 (https://doi.org/10.1038/s41467-018-04145-4)
usix = ('Lsm8', # Lsm8 is nearly identical to Lsm1, so may be some ambiguity based on recovered peptides
       'Lsm2',
       'Lsm3',
       'Lsm6',
       'Lsm5',
       'Lsm7',
       'Lsm4',
       'Prp24', # contains Rrm1-4
       'Pat1',
       'Prpf8',
       'Prpf3',
       'Prpf6',
       'Prpf4',
       'Rbm22')

# U1, Kondo et al, eLife 2015 (https://doi.org/10.7554/eLife.04986)
uone = ('Snrpc',
       'Snrpa',
       'Snrnp70',
       'Snrpb',
       'Snrpd1',
       'Snrpd2',
       'Snrpd3'
       'Snrpg',
       'Snrpe',
       'Snrpf',
       'Snrpg')

# Gene ontology splicing terms
gosplice = pd.read_csv('GO_term_splicing.csv',usecols=['Symbol'])
gosplice = list(gosplice['Symbol'])

ms_canon_rnps = (xist, 
              sevensl, 
              rmrp, 
              useven, 
              sevensk, 
              utwo, 
              usix, 
              uone,
              gosplice)

hu_canon_rnps = []

for i in ms_canon_rnps:
    hu_canon_rnps.append([u.upper() for u in i])

mousernps = dict(zip(targets, ms_canon_rnps))
humanrnps = dict(zip(targets, hu_canon_rnps))
