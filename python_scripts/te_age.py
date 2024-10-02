import os, shutil
import re
import subprocess
from Bio import SeqIO

def _get_fasta(seq_id, start, end, fasta, name):
    with open("seq.bed","w") as out:
        out.write(f"{seq_id}\t{start}\t{end}\n")
    cmd = f"bedtools getfasta -fi {fasta} -bed seq.bed > {name}"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    os.remove("seq.bed")

def get_te_age(te_index, seq_id, ltr1_coord_l, ltr2_coord_l, fasta):
    """
    LTRs coord lists with '0' based start position
    """
    if not os.path.exists(f"{te_index}_LTR_analysis/"):
        os.mkdir(f"{te_index}_LTR_analysis/")
    pwd_path = os.getcwd()
    os.chdir(f"{te_index}_LTR_analysis/")
    # get fastas
    _get_fasta(seq_id, ltr1_coord_l[0], ltr1_coord_l[1], fasta, "5LTR.fa")
    _get_fasta(seq_id, ltr2_coord_l[0], ltr2_coord_l[1], fasta, "3LTR.fa")
    # run
    ltr_len = []
    ltr_fa_list = ['5LTR.fa', '3LTR.fa']
    for ltr_fa in  ltr_fa_list:
        for r in SeqIO.parse(ltr_fa, "fasta"):
            ltr_len.append(len(r.seq))

    files2remove = ltr_fa_list
    # run stretcher
    cmd = "stretcher " + f"{ltr_fa_list[0]} {ltr_fa_list[1]} -outfile ltrIdent.txt"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # Wait for the process to finish
    process.wait()
    # parse ltr identity
    identReg = r' \(\s?(\d+\.\d+)%'
    ident = 0
    avgLtrLen = (ltr_len[0] + ltr_len[1]) / 2
    with open("ltrIdent.txt") as identFile:
        for l in identFile:
            if "# Identity:" in l:
                ident = float(re.search(identReg, l.rstrip()).group(1))
    files2remove.append("ltrIdent.txt")

    # get Kimura distance:
    ## run clustalw to get .phy file
    if not os.path.exists("outfile"):
        with open("outfile", "w") as of:
            of.write("")
    if ltr_len[0] > 0 and ltr_len[1] > 0:
        # merge two fa files:
        seq_List = []
        for r in SeqIO.parse("5LTR.fa","fasta"):
            seq_List.append(str(r.seq))
        for s in SeqIO.parse("3LTR.fa","fasta"):
            seq_List.append(str(s.seq))
        with open("infile.fa", "w") as infOut:
            infOut.write(">5LTR\n")
            infOut.write(seq_List[0] + "\n")
            infOut.write(">3LTR\n")
            infOut.write(seq_List[1] + "\n")                
        cmd1 = f"clustalw infile.fa -output=PHYLIP"
        process = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE)
        process.wait()
        
        os.rename(f'infile.phy', 'infile')
        with open("opt_dnadist", "w") as out:
            out.write("infile\nR\nD\nY\n")

        cmd1 = f"cat opt_dnadist | phylip dnadist"
        process = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE)
        process.wait()
        
        k80 = ""
        with open("outfile") as outfile:
            next(outfile)
            for l in outfile:
                lList = l.split()
                k80 = float(lList[2])  # return K80
                break
        files2remove.append("opt_dnadist")
        files2remove.append("infile")
        files2remove.append("infile.fa")
        files2remove.append("infile.dnd")
        files2remove.append("outfile")
        files2remove.append("infile")
        for file in files2remove:
            if os.path.exists(file):
                os.remove(file)
        os.chdir(pwd_path)
        shutil.rmtree(f"{te_index}_LTR_analysis/")
        return avgLtrLen, ltr_len, ident, k80