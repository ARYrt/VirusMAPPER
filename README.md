## VirusMAPPER
<!--Hi-->
![Running](/images/VirusMAPPER.gif)

Degeneration of planting material due to the accumulation of viral pathogens is a problem that affects potato crops in Colombia. In order to offer alternatives for the identification of these infectious agents, an analysis platform was developed in Python, which allows to process next generation sequence data through two alignment tools, Magic-BLAST and BLAST, also including a decision algorithm that filters the results and builds a consensus with the partial and/or complete genomes of the viral pathogens that are listed as a possible health risk according to the country's authorities.

### Supported files

This program support single or paired-end reads from Illumina sequencing technology in FASTQ format (or fastq.gz):

>@seq_ID                              [1]
>
>CTCAGCTAAATACTTTGACACCNGTANNANNNN    [2]
> 
> \+                                  [3]
>
>BBDEBDDDDHHHHFHEEEEEEEE#3AC#####   [4]

FASTQ file is composed of four lines. [1] This line includes a unique ID and information such as flow cell lane information. [2] This line includes the nucleotides of the sequence. [3] A separator item (+). [4] This line includes quality values about sequences; this quality is denominated Phred Score and is described in ASCII symbols, every simbol has a respective number asociated and a higher number signifies higher acurracy of each nucleotide.


### Pipeline

![Pipeline](/images/flow_chart2.png)


### Performance

The test data (*testfile_1.fastq* - *testfile_2.fastq*) was tested on the following computer equipment:

| Operating system | Processor| RAM | Total execution time (min) |  
| :---: | :---: | :---: | :---: |
| Windows 7 64-bit | Intel Core i5-2300 CPU @ 2,80 GHz | 8 GB | 9,71 |
| Windows 10 64-bit | Intel Core i5-6200U CPU @ 2,40 GHz | 4 GB | 18,77 |
| Windows 10 64-bit | Intel Core i5-8265U CPU @ 1,80 GHz | 8 GB | 14,24 |

----

### Installation

Sequence alignment is done through Magic-BLAST and BLAST+, the algorithms for automatic execution and filtering are written in Python.

<a name="Prerequisites"></a>
### Prerequisites

**1. BLAST+**

Please go to official [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) page and install the correct version for you operative system. To verify that the tool works correctly, enter to the terminal or the command prompt and type: 

```sh
blastn
```

The installation was **successful** if you get the following message:

```markdown
BLAST query/options error : Either a BLAST database or subject sequence(s) must be specified
Please refer to the BLAST+ user manual.
```

**2. Magic-BLAST**

Please go to official [Magic-BLAST](https://ncbi.github.io/magicblast/doc/download.html) page and install the correct version for you operative system. To verify that the tool works correctly, enter to the terminal or the command prompt and type: 

```sh
magicblast
```

The installation was **successful** if you get the following message:

```markdown
BLAST query/options error : Either a BLAST database or subject sequence(s) must be specified
Please refer to the BLAST+ user manual.
```

**3. Python**

Please go to official [Python](https://www.python.org/downloads/) page and install the correct version (3.7+) for you operative system. To verify that the tool works correctly, enter to the terminal or the command prompt and type: 

```sh
python  --version
```

The installation was **successful** if you get a message like:

```markdown
Python 3.7 (or superior)
```

<ins>Important</ins>: *If you have a version lower than Python 3, the tool may not work correctly.*

If you get an **error** like:

```markdown
"python" is not recognized as an internal or external command, operable program or batch file.
```

Use the same command, but instead of **python** try **py**. A third option is to try **python3**. If none of the previous solutions work, check the Issues section.

<ins>**More important**</ins>:

The following Python libraries are **absolutely necessary** and we recommend installing them through **pip**:

> [Pandas](https://pypi.org/project/pandas/)
> [Pillow](https://pypi.org/project/Pillow/)
> [Seaborn](https://pypi.org/project/seaborn/)
> [Matplotlib](https://pypi.org/project/matplotlib/)


Python 3.4 and later include pip by default. If it is your case enter to the terminal or the command prompt and type:

```sh
python -m pip install pandas
```
```sh
python -m pip install Pillow
```
```sh
python -m pip install seaborn
```
```sh
python -m pip install matplotlib
```

----

### Download

Before downloading the files make sure you meet the essential [prerequisites](#Prerequisites), i.e. **Python**, **BLAST+** and **Magic-BLAST**. 

You can obtain VirusMAPPER at the following [link](https://drive.google.com/drive/folders/1NBjtZpNYpZfPIJHLkxxLK43vpdDl5YPu?usp=sharing).

Once the download is complete, locate the file in your preferred folder and unzip it; when you do this you will find the following structure:

```markdown
----->VirusMAPPER.py
       
----->Databases (7 files)
       >database_viral.tsv
       >viral_refseqs.(nhr, nin, nog, nsd, nsi, nsq).
       
----->testfile_1.fastq

----->testfile_2.fastq
```

Now you can continue with the next section.

***The last update of the viral database was made on June, 2020.***

----

### Tutorial

Open the terminal or the command prompt **as administrator** and go to the folder where the downloaded files are located.

<ins>Help:</ins> If you want to know how to move between folders in the terminal or command prompt, check the following [link](https://biotecnologiamicrobianaunalmed.github.io/terminal-basics/).

You can access the program by writing the following line:

```markdown
python VirusMAPPER.py
```

Wait a moment while the user interface loads.

If you get an **error** like:

```markdown
"python" is not recognized as an internal or external command, operable program or batch file.
```

Use the same command, but instead of **python** try **py**. A third option is to try **python3**. If none of the previous solutions work, check the Issues section.

----

### Issues

- **Anaconda user's**:

> python **AND** py **AND** python3 ... is not recognized as an internal or external command, operable program or batch file.

You may have python installed through Anaconda, in this case you must open the Anaconda prompt.

### References

Boratyn GM, Thierry-Mieg J, Thierry-Mieg D, Busby B, Madden T.L. (2019) "Magic-BLAST, an accurate RNA-seq aligner for long and short reads." BMC Bioinformatics. 2019 Jul 25;20(1):405.

Zhang Z., Schwartz S., Wagner L., & Miller W. (2000), "A greedy algorithm for aligning DNA sequences" J Comput Biol 2000; 7(1-2):203-14.

### Support

Having troubles? Please contact us.
