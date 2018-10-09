![CovPlotter](http://163.172.45.124/uploads/logo.png?)


### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Plot read depth per transcript


##

![CovPlotter example](http://163.172.45.124/uploads/CovPlotter_example.png)

##

### Requirements
- python 2.7<br/>
- [yaspin](https://github.com/pavdmyt/yaspin)<br/>
- [gffutils](https://github.com/daler/gffutils)<br/>
- [pybedtools](https://github.com/daler/pybedtools)<br/>
- [samtools](https://github.com/samtools/samtools)<br/>
- tqdm

### Running
```markdown
   __      __              
  /   _   |__)| _ |_|_ _ _ 
  \__(_)\/|   |(_)|_|_(-|  
  ____________________________________________________________________

  USAGE: CovPlotter.py [OPTIONS] -g genes.txt -i listBam.txt -o output
         CovPlotter.py [OPTIONS] -n ACAD,OPA1 -i listBam.txt -o output
         CovPlotter.py [OPTIONS] -l genes.bed -i listBam.txt -o output

  [OPTIONS]
      -ref       STR      Reference hg19/hg20      (default: hg19)
      -db        STR      Database refseq/ensembl  (default: refseq)
      -tmp       DIR      Temporary folder         (default: /tmp)
      -nt        INT      Threads number           (default: 1)
      -color     BOOL     Terminal color           (true or false)
      -h,--help
      -v,--version
```