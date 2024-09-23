# Based on Neumann et al. (2019); Figure 4
# https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-018-0144-1
    
fam_pd_full_d = {"Ty3/gypsy":{"Chlamyvir":["GAG","PROT","RT","RH","INT","CHD"],
"Tcn1":["GAG","PROT","RT","RH","INT","CHD"],
"Galadriel":["GAG","PROT","RT","RH","INT","CHD"],
"Tekay":["GAG","PROT","RT","RH","INT","CHD"],
"Reina":["GAG","PROT","RT","RH","INT","CHD"],
"CRM":["GAG","PROT","RT","RH","INT","CHDCR"],
"Phygy":["GAG","PROT","RT","RH","INT"],
"Selgy":["GAG","PROT","RT","RH","INT"],
"Athila":["GAG","PROT","RT","RH","INT"],
"TatI":["GAG","PROT","aRH""RT","RH","INT"],
"TatII":["GAG","PROT","aRH""RT","RH","INT"],
"TatIII":["GAG","PROT","RT","RH","INT","aRH"],
"Ogre":["GAG","PROT","RT","RH","aRH","INT"],
"Retand":["GAG","PROT","RT","RH","aRH","INT"]},
"Ty1/copia":{"Osser":["GAG","PROT","INT","RT","RH"],
"Bryco":["GAG","PROT","INT","RT","RH"],
"Lyco":["GAG","PROT","INT","RT","RH"],
"Gymco-I":["GAG","PROT","INT","RT","RH"],
"Gymco-II":["GAG","PROT","INT","RT","RH"],
"Gymco-III":["GAG","PROT","INT","RT","RH"],
"Gymco-IV":["GAG","PROT","INT","RT","RH"],
"Ale":["GAG","PROT","INT","RT","RH"],
"Ivana":["GAG","PROT","INT","RT","RH"],
"Ikeros":["GAG","PROT","INT","RT","RH"],
"Tork":["GAG","PROT","INT","RT","RH"],
"Alesia":["GAG","PROT","INT","RT","RH"],
"Angela":["GAG","PROT","INT","RT","RH"],
"Bianca":["GAG","PROT","INT","RT","RH"],
"SIRE":["GAG","PROT","INT","RT","RH"],
"TAR":["GAG","PROT","INT","RT","RH"]}}

def eval_te_autonomy(sfam, fam, pd_list):
    out_pd_list = []
    for dom in fam_pd_full_d[sfam][fam]:
        if dom in pd_list:
            out_pd_list.append(dom)
        else:
            out_pd_list.append("NA")
    autonomy_status = "autonomous"
    if "NA" in out_pd_list:
        autonomy_status = "non_autonomous"
    return autonomy_status, "|".join(out_pd_list)