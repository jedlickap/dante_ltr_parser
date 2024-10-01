# Based on Neumann et al. (2019); Figure 4
# https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-018-0144-1
    
fam_pd_full_d = {'Ty1/copia': {('GAG', 'PROT', 'INT', 'RT', 'RH'): ['Osser',
                                                    'Bryco',
                                                    'Lyco',
                                                    'Gymco-I',
                                                    'Gymco-II',
                                                    'Gymco-III',
                                                    'Gymco-IV',
                                                    'Ale',
                                                    'Ivana',
                                                    'Ikeros',
                                                    'Tork',
                                                    'Alesia',
                                                    'Angela',
                                                    'Bianca',
                                                    'SIRE',
                                                    'TAR']},
 'Ty3/gypsy': {('GAG', 'PROT', 'RT', 'RH', 'INT'): ['Phygy', 'Selgy', 'Athila'],
               ('GAG', 'PROT', 'RT', 'RH', 'INT', 'CHD'): ['Chlamyvir',
                                                           'Tcn1',
                                                           'Galadriel',
                                                           'Tekay',
                                                           'Reina'],
               ('GAG', 'PROT', 'RT', 'RH', 'INT', 'CHDCR'): ['CRM'],
               ('GAG', 'PROT', 'RT', 'RH', 'INT', 'aRH'): ['TatIII'],
               ('GAG', 'PROT', 'RT', 'RH', 'aRH', 'INT'): ['Ogre', 'Retand'],
               ('GAG', 'PROT', 'aRH','RT', 'RH', 'INT'): ['TatI', 'TatII']}}


def eval_te_autonomy(sfam, fam, pd_list):
    out_pd_list = []
    for dom in [dom_pat for dom_pat in fam_pd_full_d[sfam] if fam in fam_pd_full_d[sfam][dom_pat]][0]:
        if dom in pd_list:
            out_pd_list.append(dom)
        else:
            out_pd_list.append("NA")
    autonomy_status = "autonomous"
    if "NA" in out_pd_list:
        autonomy_status = "non_autonomous"
    return autonomy_status, "|".join(out_pd_list)