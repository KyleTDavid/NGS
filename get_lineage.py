from ete3 import NCBITaxa
import sys
ncbi = NCBITaxa()

d={}

print 'taxid\tSpecies\tGenus\tFamily\tOrder\tClass\tPhylum\tKingdom'


with open(sys.argv[1]) as f:
    for line in f:
        line=line.lstrip().rstrip()
        ID=line
        Species='NA'
        Genus='NA'
        Family='NA'
        Order='NA'
        Class='NA'
        Phylum='NA'
        Kingdom='NA'
        if line == 'NA':
        	continue
        else:
			for clade in ncbi.get_lineage(line):
				if ncbi.get_rank([clade]).values()[0]=='species':
					Species=ncbi.get_taxid_translator([clade]).values()[0]
				if ncbi.get_rank([clade]).values()[0]=='genus':
					Genus=ncbi.get_taxid_translator([clade]).values()[0]
				if ncbi.get_rank([clade]).values()[0]=='family':
					Family=ncbi.get_taxid_translator([clade]).values()[0]
				if ncbi.get_rank([clade]).values()[0]=='order':
					Order=ncbi.get_taxid_translator([clade]).values()[0]
				if ncbi.get_rank([clade]).values()[0]=='class':
					Class=ncbi.get_taxid_translator([clade]).values()[0]
				if ncbi.get_rank([clade]).values()[0]=='phylum':
					Phylum=ncbi.get_taxid_translator([clade]).values()[0]
				if ncbi.get_rank([clade]).values()[0]=='kingdom':
					Kingdom=ncbi.get_taxid_translator([clade]).values()[0]
        print ID,'\t',Species,'\t',Genus,'\t',Family,'\t',Order,'\t',Class,'\t',Phylum,'\t',Kingdom




