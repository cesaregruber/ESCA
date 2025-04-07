

import optparse
from Bio import SeqIO
# import IUPACfunc




#
#from os import listdir, chdir
#from os.path import isfile, isdir, join
#from shutil import copyfile
#import subprocess



def translate(seq): 
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein+= table[codon] 
    return protein 





def variantdefine(varlist, infasta,  annotazioni):
    #infasta="/usr/bin/SCV2AssemblyPackage/NC_045512.2.fasta" #/home/cesare/SARSCoV2/AnalisiClustersPisanaLatina/settembre/
    #annotazioni="/usr/bin/SCV2AssemblyPackage/NC_045512.2 AnnotationsMOD.tsv.csv" #/home/cesare/SARSCoV2/AnalisiClustersPisanaLatina/settembre/
    seqdic= SeqIO.to_dict(SeqIO.parse(infasta, "fasta"))
    refseq=str(seqdic['NC_045512.2'].seq.upper())
    regions={}
    indata=open(annotazioni, 'r')
    casa = 0
    for line in indata :
        ls = line.split("\t")
        if "gene" in line :
            name=ls[0].split('gene')[0][:-1]
            startpos = int(ls[2])
            endpos = int(ls[3])
            if (endpos - startpos + 1) != int(ls[4]) :
                casa = 2
#                print "Errore nel gene: " 
#                print name
#                print ls[3],  ls[2],  ls[4]
            myNNseq =   refseq[startpos-1:endpos]
            myAAseq = translate(myNNseq)
            regions[name] = [startpos, endpos, myNNseq, myAAseq]
    
    variants=[]
    if isinstance(varlist, str) :
        variants=varlist.split(',')
    elif isinstance(varlist, list) : 
        variants = varlist
    else :
        casa = 3
#        print "ERORRE LISTA VARIANTI NON RICONOSCIUTA!"
    
    AAvariants = []
    
    for v in variants :
        prevN = v[0] 
        nextN = v[-1]
        posN = int(v[1:-1])
        if refseq[posN-1] != prevN :
            casa = 4
 #           print "ERRORE: variante " + v + " non trovata!"
        else :
            newseq = refseq[:posN-1] + nextN + refseq[posN:] 
        
        trovata = 0
            
        for gene in regions :
            if (posN  > regions[gene][0] ) and ( posN  < regions[gene][1]  ) :
                trovata = 1
                casa = 5
#                print "La varinte " + v + " si trova nel gene " + gene 
                sp = regions[gene][0] - 1
                ep = regions[gene][1]
                newseqgeneN=newseq[sp:ep]
                oldseqgeneN = regions[gene][2] 
                if newseqgeneN == oldseqgeneN :
                    casa = 1
#                    print "Errore: variante " + v + " non trovata  nel gene " + gene 
                else :
                    newseqgeneA = translate(newseqgeneN)
                    oldseqgeneA = regions[gene][3] 
                    if newseqgeneA == oldseqgeneA :
#                        print "La variante " + v + " e sinonima"
                        AAvar = gene + ':'+ "Synonymous"
                    else :
                        for a in range(len(newseqgeneA)) :
                            if newseqgeneA[a] != oldseqgeneA[a] :
                                AAvar = gene + ':'+ oldseqgeneA[a] + str(a+1) + newseqgeneA[a] 
#                                print "Ecco la variante :"
#                                print AAvar
                                
                    AAvariants.append(AAvar)
        
        if trovata == 0 :
#            print     "La variante " + v + " non si trova in una regione codificante"
            AAvar = 'NC:NC'
            AAvariants.append(AAvar)
            
    return(AAvariants)
        


# funzione che consente di trovare le varianti 
#  fra seq1 e seq2 , registrando le posizioni di seq1
def trovavarianti(seq1, seq2):
    sequ1=seq1.upper()
    sequ2=seq2.upper()
    caratteracci=set( ['-','M','R','W','S','Y','K','V','H','D','B','N', '~', '.', ','])
    varianti = {}
    
    if len(sequ1) == len(sequ2) :
        position = 0
        
        for pair in zip(sequ1,sequ2):
            position +=1 
            if (pair[0] != pair[1]) and ( caratteracci & set(pair)  == set() )  :          #controllo che la variante sia fra due nucleoditi reali
                if pair[1] in 'ACGT' :
                    varianti[position] = pair
                
        return varianti





# INIZIO MAIN !!!! ####################################################################################################


#
#infasta="/home/cesare/SARSCoV2/LavoroBangladeshMartina/NuoveSARScov2consensus_200715_SecondoChip_SoloBangla.afa"#"NuoveSARScov2consensusUpdate_200828.afa"
##
#outvariants=open("/home/cesare/SARSCoV2/LavoroBangladeshMartina/NuoveSARScov2consensus_200715_SecondoChip_SoloBangla_varianti.csv", 'w') #"NuoveSARScov2consensusUpdate_200828noIO.csv"
#
#



parser = optparse.OptionParser(description="QUESTO PROGRAMMA serve a trovare tutte le varianti fra le nuove sequenze ed il riferimento Wuhan-Hu1")
parser.add_option('-i', dest='in_fasta', help='file fasta allineato Con la Sequenza Wuhan-Hu-1 versione NC_045512 di NCBI di input ')
parser.add_option('-o', dest='out_csv', help='nome csv output')
parser.add_option('-r', dest='in_ref', help='nome file referenza')
parser.add_option('-a', dest='in_ann', help='nome file annotazione referenza')

args, options = parser.parse_args()

infasta=open(args.in_fasta,'r')
outvariants=open(args.out_csv,'w')
refasta=str(args.in_ref)
annotazioni=str(args.in_ann)



seqs= SeqIO.to_dict(SeqIO.parse(infasta, "fasta"))

names = seqs.keys()
names.sort()

initseq=[0]*(len(names) -1)

#newnames=[k.split('_')[0] for k in names]


myref = seqs['NC_045512.2']

myrefseq = str(myref.seq).upper().replace('U','T')

Allvariants = {}

variantpositions = set([])
nuclvars= set([])

for name in names :
    if name != myref.id : 
        namepos =  names.index(name) - 1 
        nameseq = str(seqs[name].seq).upper()
        myvariants = trovavarianti( myrefseq, nameseq )
        
        variantpositions =  variantpositions | set(myvariants.keys())
        
        mynuclvars = []
        for var in myvariants :
           myvar = myvariants[var][0]+str(var)+myvariants[var][1]
           mynuclvars.append(myvar)
           
        nuclvars = nuclvars | set(mynuclvars)
        
        myname = name.split('.')[0] 
        
        Allvariants[myname] = myvariants
        
        
        
        

ordnuclvars = sorted(list(nuclvars),key=lambda x:int(x[1:-1]) )

aminoacidvars = variantdefine(ordnuclvars, refasta,  annotazioni)

#print aminoacidvars



mynames = [name.split('.')[0] for name in names]
mynames.remove("NC_045512")

variantposlist= list(variantpositions)
variantposlist.sort()

# if len(aminoacidvars) != variantposlist :
    # print "ERRORE: numero di varianti nucleotidi diverso dal numero di varianti analizzate per la traduzione"
    # print variantposlist
    # print aminoacidvars
    
    
outvariants.write("POSITION,Wuhan-Hu-1,AminoVar,"+','.join( mynames ) + '\n'  )

i=0
for pos in  variantposlist :
    
    refbase= myrefseq[pos-1].upper()
    myAAvar= aminoacidvars[i]
    
    i=i+1
    outvariants.write(str(pos)+','+refbase + ','+myAAvar+',')
    
    for myname in mynames :
        if pos not in Allvariants[myname]  :
            trovata = 0
            for name in names : 
                if myname in name :
                    nameseq = str(seqs[name].seq).upper()
                    if nameseq[pos-1] != refbase : 
                        outvariants.write( nameseq[pos-1] + ',')
                        trovata = 1
                    else :
                        outvariants.write(',')
                        trovata = 1
            if trovata == 0 :
                print "ERRORE!!! Ambiguita non risolta: ",   myname ,  str(pos)
        else :
            controlbase = Allvariants[myname][pos][0] 
            if controlbase != refbase :
                print 'ERROREEEEEE !!!!!!!',  myname ,  str(pos)
            else  :
                outvariants.write( Allvariants[myname][pos][1] + ',')

    outvariants.write('\n')
    
    
outvariants.close()









