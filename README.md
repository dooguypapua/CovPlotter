![CovPlotter](http://163.172.45.124/uploads/logo.png?)


### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Plot read depth per transcript


##

![CovPlotter example](http://163.172.45.124/uploads/CovPlotter_example.png)

##


### Requirements
- python3<br/>
- [wkhtmltopdf](https://github.com/wkhtmltopdf/wkhtmltopdf)<br/>
- python package: [Pillow](https://github.com/python-pillow/Pillow), [gffutils](https://github.com/daler/gffutils), [pybedtools](https://github.com/daler/pybedtools), [binaryornot](https://github.com/audreyr/binaryornot), [yaspin (0.12.0)](https://github.com/pavdmyt/yaspin), [tqdm](https://github.com/tqdm/tqdm)


### Installation
```markdown
wget https://github.com/wkhtmltopdf/wkhtmltopdf/releases/download/0.12.4/wkhtmltox-0.12.4_linux-generic-amd64.tar.xz
tar -xvf wkhtmltox-0.12.4_linux-generic-amd64.tar.xz
cd wkhtmltox
sudo cp -R ./* /usr/
git clone https://github.com/dooguypapua/CovPlotter.git
cd CovPlotter
pip3 install -r requirements.txt
```

### Usage
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