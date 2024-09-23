import os
import re
import subprocess
from Bio import SeqIO

def _get_fasta(seq_id, start, end, fasta, name):
    with open("seq.bed","w") as out:
        out.write(f"{seq_id}\t{start}\t{end}\n")
    cmd = f"bedtools getfasta -fi {fasta} -bed seq.bed > {name}"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    print(f"BEDTOOLS GETFASTA process.returncode: {process.returncode}")
    os.remove("seq.bed")

def get_te_age(seq_id, ltr1_coord_l, ltr2_coord_l, fasta):
    """
    LTRs coord lists with '0' based start position
    """
    # get fastas
    _get_fasta(seq_id, ltr1_coord_l[0], ltr1_coord_l[1], fasta, "5LTR.fa")
    _get_fasta(seq_id, ltr2_coord_l[0], ltr2_coord_l[1], fasta, "3LTR.fa")

    # run
    ltr_len = []
    for ltr_fa in  ['5LTR.fa','3LTR.fa']:
        for r in SeqIO.parse(ltr_fa, "fasta"):
            ltr_len.append(len(r.seq))

    files2remove = ['5LTR.fa','3LTR.fa']
    # run stretcher
    cmd = "stretcher " + "5LTR.fa 3LTR.fa -outfile ltrIdent.txt"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    print(f"STRETCHER process.returncode: {process.returncode}")
    # parse ltr identity
    identReg = r' \(\s?(\d+\.\d+)%'
    ident = 0
    avgLtrLen = (ltr_len[0] + ltr_len[1]) / 2
    with open("ltrIdent.txt") as identFile:
        for l in identFile:
            if "# Identity:" in l:
                print(l.rstrip())
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
        cmd1 = "clustalw infile.fa -output=PHYLIP"
        process = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(f"CLUSTALW with phylip output process.returncode: {process.returncode}")

        os.rename('infile.phy', 'infile')
        with open("opt_dnadist", "w") as out:
            out.write("infile\nR\nD\nY\n")

        cmd1 = "cat opt_dnadist | phylip dnadist"
        process = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(f"PHYLIP K80 divergence process.returncode: {process.returncode}")

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
        return avgLtrLen, ltr_len, ident, k80