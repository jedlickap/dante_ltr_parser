import pandas as pd
import subprocess

style = """
<style>
body {
  font-family: Arial, sans-serif; /* Set font to Arial for the entire document */
}
table, th, td {
  border: 1px solid black;
  border-collapse: collapse;
}

th, td {
  padding: 5px;
}
</style>
"""

def generate_report(out_path, species):
    out_folder_path = out_path + "/tabs_and_plots/"
    df = pd.read_csv(f"{out_folder_path}summary_table.csv")

    # Apply the style to the rows where 'te_fam' contains 'Summary'
    def highlight_summary(row):
        if "Summary" in str(row['te_fam']):
            return ['background-color: lightgrey'] * len(row)
        else:
            return [''] * len(row)

    # Use the Styler to apply the highlighting
    styled_df = (df.style
                 .apply(highlight_summary, axis=1)
                 .format(precision=2))  # Set float precision to 2 decimals

    # Convert to HTML with the applied style
    sum_tab = styled_df.hide(axis="index").to_html()

    # Save or return the HTML table (saving in this example)
    with open(f"{out_folder_path}styled_summary_table.html", "w") as f:
        f.write(sum_tab)

    md = f"""

# LTR retrotransposons in *{species}*

### [visualization of DANTE_LTR GFF3 output]

## Contents

  - [Summary table](#summary-table)
  - [Quality of detected TEs](#quality-of-detected-tes)
  - [Proportion of autonomous and non-autonomous TEs](#proportion-of-autonomous-and-non-autonomous-tes)
  - [TEs and LTR lengths](#tes-and-ltr-lengths)
  - [LTR identity and TE age](#ltr-identity-and-te-age)
  - [TE distribution in chromosomes](#te-distribution-in-chromosomes)
  - [TE age distribution in chromosomes](#te-age-distribution-in-chromosomes)
  - [TE family specific age distribution in chromosomes](#te-family-specific-age-distribution-in-chromosomes)

## Summary table

{style}

{sum_tab}

## Quality of detected TEs

<img src="./tabs_and_plots/tsd_pbs_summary.png" alt="TEs quality" height="500" width="900"/>

## Proportion of autonomous and non-autonomous TEs

<img src="./tabs_and_plots/autnom_nonauton_summary.png" alt="TE autonomy" height="500" width="900"/>

## TEs and LTR lengths

<img src="./tabs_and_plots/te_ltr_lengths.png" alt="TE & LTR lengths" height="500" width="900"/>

## LTR identity and TE age

<img src="./tabs_and_plots/ltr_ident_mya.png" alt="LTR identity & TE age" height="500" width="900"/>

## TE distribution in chromosomes

<img src="./tabs_and_plots/te_cnt_chr_specific_density.png" alt="TEs in chromosomes" height="1200" width="1200"/>

## TE age distribution in chromosomes

<img src="./tabs_and_plots/te_age_chr_specific_density.png" alt="TE age in chromosomes" height="1200" width="1200"/>

## TE family specific age distribution in chromosomes

<img src="./tabs_and_plots/te_age_chr_fam_spec_density.png" alt="TE age family-specific distribution" height="1900" width="1900"/>

    """
    with open(f"{out_path}/report.md", "w") as textfile:
        textfile.write(md)

    # md to html conversion in pandoc
    # pandoc -f markdown report.md > report_pandoc.html
    cmd = f"pandoc -f markdown {out_path}/report.md > {out_path}/report.html"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    print(f"PANDOC MD2HTML CONVERSION process.returncode: {process.returncode}")
