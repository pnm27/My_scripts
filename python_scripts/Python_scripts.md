**Script Description**
<h2>table_browser_to_SCRC.py</h2>

This script accepts text files containing fasta records (see [here]() for toy input) to produce a `txt` file that can serve as input to [discovery]() mode of [SCReadCounts](https://github.com/HorvathLab/NGS/tree/master/SCReadCounts)

This script reads the input fasta records and for each "reference" nucleotide at each position of the sequence it uses the remaining 3 nucleotides as the "alternate". Check the [toy]() inputs for more description of the script.