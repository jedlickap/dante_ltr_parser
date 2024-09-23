from python_scripts.te_age import get_te_age
from python_scripts.te_position_in_bins import gen_all_bin_lists, position_in_bins
from python_scripts.eval_te_autonomy import eval_te_autonomy
from python_scripts.GFFsParsers import DanteLTrGFF
from pathlib import Path

class TEAnalyzer:
    def __init__(self, dnt_gff, genome_fa, ssr, out_path):
        self.dnt_gff = dnt_gff
        self.genome_fa = genome_fa
        self.ssr = ssr
        self.out_path = out_path
        self.dntObj = DanteLTrGFF(dnt_gff)
        self.dnt_dict = self.dntObj.process_data()
        self.fasta_bin_dict = gen_all_bin_lists(genome_fa)

    def _get_te_fam(self, final_classif):
        ll = final_classif.split("|")
        for lab in ll:
            if lab.startswith("Ty"):
                sfam = lab
        fam = ll[-1]
        return sfam, fam

    def _get_prot_doms(self, prot_dom_list):
        return [pd.attributes['Name'] for pd in prot_dom_list]

    def _get_age_cat(self, mya):
        age_dict = {"A": [0, 1], "B": [1, 2], "C": [2, 3], "D": [3, 4], "E": [4, 5], "F": [5, 10], "G": [10, 1000]}
        for cat, (s, e) in age_dict.items():
            if s <= mya < e:
                return cat

    def csv_generator(self):
        out_csv_path = f"{self.out_path}/{Path(self.dnt_gff).stem}_TE_characteristics.csv"
        with open(out_csv_path, "w") as csv:
            csv.write("chromosome,te_id,te_sfam,te_fam,tsd,pbs,prot_doms,autonomy_stat,te_length,ltr_avg_len,ltr5_len,ltr3_len,ltr_identity,K80,MYA,age_cat,bin_10kbp,bin_100kbp,bin_1Mbp,inters_10kbp,inters_100kbp,inters_1Mbp\n")
            for te in self.dnt_dict:
                seq_id = self.dnt_dict[te]['transposable_element'][0].seqid
                te_id = self.dnt_dict[te]['transposable_element'][0].attributes['ID']
                if "partial" not in te_id:
                    print(f"CHROMOSOME: {seq_id}\nTEID: {te_id}")
                    te_start = self.dnt_dict[te]['transposable_element'][0].start
                    te_end = self.dnt_dict[te]['transposable_element'][0].end
                    te_len = (te_end - te_start) + 1
                    tsd = bool(self.dnt_dict[te]['target_site_duplication'])
                    pbs = bool(self.dnt_dict[te]['primer_binding_site'])
                    ltr1_coord_l = [self.dnt_dict[te]['long_terminal_repeat'][0].start, self.dnt_dict[te]['long_terminal_repeat'][0].end]
                    ltr2_coord_l = [self.dnt_dict[te]['long_terminal_repeat'][1].start, self.dnt_dict[te]['long_terminal_repeat'][1].end]
                    avgLtrLen, ltr_len, ident, k80 = get_te_age(seq_id, ltr1_coord_l, ltr2_coord_l, self.genome_fa)
                    ltr_len_str = ",".join(map(str, ltr_len))
                    mya = round((k80 / (2 * self.ssr)) / 1000000, 4) if k80 else "NA"
                    age_cat = self._get_age_cat(mya) if mya != "NA" else "NA"
                    sfam, fam = self._get_te_fam(self.dnt_dict[te]['transposable_element'][0].attributes['Final_Classification'])
                    pd_list = self._get_prot_doms(self.dnt_dict[te]['protein_domain'])
                    autonomy_status, prot_doms = eval_te_autonomy(sfam, fam, pd_list)
                    bin_d = position_in_bins(seq_id, te_start, te_end, te_len, self.fasta_bin_dict)
                    bin_sec_string = ",".join([bin_d[k][0].split("|")[0] for k in bin_d for rec in bin_d[k]])
                    bin_inters_string = ",".join([bin_d[k][0].split("|")[1] for k in bin_d for rec in bin_d[k]])
                    csv.write(f"{seq_id},{te_id},{sfam},{fam},{tsd},{pbs},{prot_doms},{autonomy_status},{te_len},{avgLtrLen},{ltr_len_str},{ident},{k80},{mya},{age_cat},{bin_sec_string},{bin_inters_string}\n")
                    print("\n")
        return out_csv_path

# Example usage
# if __name__ == "__main__":
#     args = args_from_parser()
#     analyzer = TEAnalyzer(args.dante_ltr_gff, args.genome_fasta, args.synonymous_substitution_rate, args.output_path)
#     out_csv_path = analyzer.csv_generator()
