from Bio import SeqIO

def get_mbp_bins(head, ch_len):
    l = []
    step = 1
    start = 0
    end = 1_000_000
    while (ch_len - end) >= 1_000_000:
        l.append(f"{head}|{step}|{start}|{end}")
        start = end
        end += 1_000_000
        step += 1
    l.append(f"{head}|{step}|{start}|{end}")
    return l

def get_100kbp_bins(head, ch_len):
    l = []
    step = 1
    start = 0
    end = 100000
    while (ch_len - end) >= 100000:
        l.append(f"{head}|{step}|{start}|{end}")
        start = end
        end += 100000
        step += 1
    l.append(f"{head}|{step}|{start}|{end}")
    return l

def get_10kbp_bins(head, ch_len):
    l = []
    step = 1
    start = 0
    end = 10000
    while (ch_len - end) >= 10000:
        l.append(f"{head}|{step}|{start}|{end}")
        start = end
        end += 10000
        step += 1
    l.append(f"{head}|{step}|{start}|{end}")
    return l

def get_fasta_dict(fasta):
    # genomic fasta statistics
    fasta_dict = {}
    for rec in SeqIO.parse(fasta,"fasta"):
        fasta_dict[rec.id] = len(rec.seq)
    return fasta_dict

def get_intersect(bin_s, bin_e, te_start, te_end):
    # Calculate the start and end points of the intersection
    start_intersect = max(bin_s, te_start)
    end_intersect = min(bin_e, te_end)
    
    # If there's no intersection, return 0
    if start_intersect >= end_intersect:
        return 0
    
    # Calculate the size of the intersection
    intersect_size = (end_intersect - start_intersect) + 1
    
    return intersect_size

def gen_all_bin_lists(fasta):
    bin_d = {"bin10k":{},"bin100k":{},"bin1M":{}}
    fasta_dict = get_fasta_dict(fasta)
    for seq_id in fasta_dict:
        bin_d["bin10k"][seq_id] = get_10kbp_bins(seq_id, fasta_dict[seq_id])
        bin_d["bin100k"][seq_id] = get_100kbp_bins(seq_id, fasta_dict[seq_id])
        bin_d["bin1M"][seq_id] = get_mbp_bins(seq_id, fasta_dict[seq_id])
    return bin_d

def position_in_bins(seq_id, te_start, te_end, te_len, fasta_bin_dict):
    bin_d = {"bin10k":[],"bin100k":[],"bin1M":[]}
    for bin_cat in fasta_bin_dict:
        for seq_head in fasta_bin_dict[bin_cat]:
            for bin in fasta_bin_dict[bin_cat][seq_head]:
                bin_l = bin.split("|")
                b_sec, bin_s, bin_e = bin_l[1], int(bin_l[2]), int(bin_l[3])
                if get_intersect(bin_s, bin_e, te_start, te_end):
                    intersect_len = get_intersect(bin_s, bin_e, te_start, te_end)
                    if not bin_d[bin_cat]: 
                        bin_d[bin_cat].append(f"{b_sec}|{intersect_len}")
    return bin_d