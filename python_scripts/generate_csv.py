from python_scripts.te_age import get_te_age
from python_scripts.te_position_in_bins import gen_all_bin_lists, position_in_bins
from python_scripts.eval_te_autonomy import eval_te_autonomy
from python_scripts.fasta_parser import FastaParser
from python_scripts.GFFsParsers import DanteLTrGFF
from pathlib import Path
import os
import multiprocessing

class TEAnalyzer:
    def __init__(self, dnt_gff, genome_fa, ssr, out_path, threads):
        self.dnt_gff = dnt_gff
        self.genome_fa = genome_fa
        self.ssr = ssr
        self.out_path = out_path
        self.threads = threads
        self.dntObj = DanteLTrGFF(dnt_gff)
        self.dnt_dict = self.dntObj.process_data()
        self.fasta_bin_dict = gen_all_bin_lists(genome_fa)
	# Get length of thte longest sequence
        self.fasta_parser = FastaParser(self.genome_fa)
        self.length_in_mb = self.fasta_parser.get_longest_sequence_in_mb()
        
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
                
    def _get_csv_line(self, te, te_index):
        seq_id = self.dnt_dict[te]['transposable_element'][0].seqid
        te_id = self.dnt_dict[te]['transposable_element'][0].attributes['ID']        
        te_start = self.dnt_dict[te]['transposable_element'][0].start
        te_end = self.dnt_dict[te]['transposable_element'][0].end
        te_len = (te_end - te_start) + 1
        tsd = bool(self.dnt_dict[te]['target_site_duplication'])
        pbs = bool(self.dnt_dict[te]['primer_binding_site'])
        ltr1_coord_l = [self.dnt_dict[te]['long_terminal_repeat'][0].start, self.dnt_dict[te]['long_terminal_repeat'][0].end]
        ltr2_coord_l = [self.dnt_dict[te]['long_terminal_repeat'][1].start, self.dnt_dict[te]['long_terminal_repeat'][1].end]
        avgLtrLen, ltr_len, ident, k80 = get_te_age(te_index, seq_id, ltr1_coord_l, ltr2_coord_l, self.genome_fa)
        ltr_len_str = ",".join(map(str, ltr_len))
        mya = round((k80 / (2 * self.ssr)) / 1000000, 4) if k80 > 0 else 0.0
        age_cat = self._get_age_cat(mya) if mya != "NA" else "NA"
        sfam, fam = self._get_te_fam(self.dnt_dict[te]['transposable_element'][0].attributes['Final_Classification'])
        pd_list = self._get_prot_doms(self.dnt_dict[te]['protein_domain'])
        autonomy_status, prot_doms = eval_te_autonomy(sfam, fam, pd_list)
        csv_lines = []
        bin_d = position_in_bins(seq_id, te_start, te_end, te_len, self.fasta_bin_dict)
        for bin_res in bin_d:
            bin_sec = bin_d[bin_res][0].split("|")[0]
            bin_inters = bin_d[bin_res][0].split("|")[1]
            csv_lines.append(f"{bin_res},{seq_id},{te_id},{sfam},{fam},{tsd},{pbs},{prot_doms},{autonomy_status},{te_len},{avgLtrLen},{ltr_len_str},{ident},{k80},{mya},{age_cat},{bin_sec},{bin_inters}\n")
        return csv_lines

    # Helper function to unpack the tuple arguments
    def _get_csv_line_helper(self, data):
        te, te_index, tes_total = data
        return self._get_csv_line(te, te_index)

    def _assign_category(self, seq_len):
        if seq_len <= 1:
            return 'bin10k'
        elif 1 < seq_len <= 100:
            return 'bin100k'
        else:
            return 'bin1M'
    
    def csv_generator(self):
        out_csv_path = f"{self.out_path}/{Path(self.dnt_gff).stem}_TE_characteristics.csv"
        
        # Filter list of TE IDs for non-partial TEs
        teid_list = [te for te in self.dnt_dict if "TE_partial" not in self.dnt_dict[te]['transposable_element'][0].attributes['ID']]

        # Prepare the data to pass to multiprocessing (tuples of (te, te_index))
        te_data = [(te, index, len(teid_list)) for index, te in enumerate(teid_list)]

        # Get number of available CPU cores and limit the number of threads (processes)
        num_workers = self.threads  # Use half the number of cores, at least 1 worker

	    # Create a pool of workers with limited threads
        with multiprocessing.Pool(processes=num_workers) as pool:
            # Use pool.map to parallelize the _get_csv_line function
            csv_lines = []
            for _ in pool.imap_unordered(self._get_csv_line_helper, te_data, chunksize=20):
            # csv_lines = pool.map(self._get_csv_line_helper, te_data)
                csv_lines.extend(_)
                # report the number of remaining tasks
                processed_tes = len({line.split(",")[1] + "|" + line.split(",")[2] for line in csv_lines})
                print(f'{(processed_tes)}/{len(teid_list)} TEs completed')

        # Write to the CSV file
        for bin_res in ['bin10k','bin100k','bin1M']:
            csv_bin_path = out_csv_path.replace(".csv",f"_{bin_res}.csv")
            write_header = not os.path.exists(csv_bin_path)
            
            with open(csv_bin_path, "a") as h:
                if write_header:
                    h.write("chromosome,te_id,te_sfam,te_fam,tsd,pbs,prot_doms,autonomy_stat,te_length,"
                        "ltr_avg_len,ltr5_len,ltr3_len,ltr_identity,K80,MYA,age_cat,"
                        f"bin_{bin_res},inters_{bin_res}\n")
                
                for line in csv_lines:
                    if line.startswith(bin_res):
                        line = ",".join(line.split(",")[1:])
                        h.write(line + "\n")

        # Return path for CSV with optimal bin category
        bin_cat = self._assign_category(self.length_in_mb)
        return out_csv_path.replace(".csv",f"_{bin_cat}.csv")
