

import IUPACfunc
import optparse

import itertools
import numpy  


#
#
#from os import chdir
#chdir("/home/cesare/SARSCoV2/NuovaAnalisiGeneSconAmpliconi/workdir_fineAprile/Pt5BAL")
#
#infile=open("Pt5BAL_QualFilTrimmed_ONLY_WuhanHu1_Sorted.mpileup",'r')
#outfile=open("Pt5BAL_QualFilTrimmed_ONLY_WuhanHu1_Sorted_TutteVariantiPerDiversitaNEWQUASPROVA.csv",'w')
#outfasta=open("Pt5BAL_consensus_NewQuasPROVA.fasta",'w')
#
#outquas=open( "Pt5BAL_quaispecies_NewQuasPROVA.fasta",'w')       
#
#fastaheader="Pt5BAL_consensus_NewQuasPROVA"
#
#outfasta.write(">"+ "Pt5BAL_consensus_NewQuasPROVA" +"\n")
#





parser = optparse.OptionParser(description="QUESTO PROGRAMMA serve a convertire un file mpileup in un file csv e ricostruire la sequenza consensus")
parser.add_option('-i', dest='in_pileup', help='file mpileup di input')
parser.add_option('-o', dest='out_csv', help='nome csv output')
parser.add_option('-f', dest='out_fa', help='nome fasta output')
#parser.add_option('-q', dest='out_qu', help='nome fasta output quasispecie')
parser.add_option('-n', dest='name', default="Consensus" ,  help='name header fasta')



args, options = parser.parse_args()


infile = open(args.in_pileup,'r')
outfile = open(args.out_csv,'w')
outfasta=open( args.out_fa,'w')
#outquas=open( args.out_qu,'w')                  

fastaheader = args.name 

outfasta.write(">"+ fastaheader + '_consensus\n')







allelecountKeys=['A', 'C', 'G', 'T', 'N', 'insertions','deletions' ]


ambidic ,  ambinonambidic ,  invambidic ,  complambidic ,  ambikeys  ,  ambiset  = IUPACfunc.iupacdic()



outfile.write('\t'.join(['position', 'Ref/consensus',  'coverage', 'A', 'C', 'G', 'T', 'N', 'insertions','deletions' ]) + '\n' )



# PARAMETRO :
minallelefreq=0.05

newseq=''
positions=[]

quasalls = []
quasfreqs = []
quaspos = []

for line in infile :
    ls = line.split('\t')
    
    position = ls[1]
    refallele = ls[2]
    coverage = ls[3]
    
    string = ls[4]
    
    

    allelecount = {}
    
    allelecount['A'] = 0
    allelecount['C'] = 0
    allelecount['G'] = 0
    allelecount['T'] = 0
    allelecount['N'] = 0
    allelecount['insertions'] = 0
    allelecount['deletions'] = 0

    
    ref = refallele
    
    altc = ''
    altf=[]
    
    newallele = ""
    indel='false'
    indelcount=0
    indelmax=0
    
    insert = ''
    
    for c in string :
        if c == '+' :
            indel = 'true'
            insert = 'true'
        
        elif c == '-' :
            indel = 'true'
            insert = 'false'

        elif c in ".," and indel=='false' :  # for il punto rev la virgola
            allelecount[ref] += 1
        
        elif c in 'AGTCNagtcn' and indel=='false'   :
            newallele = c.upper()
            allelecount[newallele] += 1
        
        elif (c in '0123456789') and (indel == 'true') and (indelcount == 0 ) :
            indelcount = int(c)
            indelmax = indelcount
        
                
        elif c == '*' :
            allelecount['deletions'] += 1
        
        elif (c in '0123456789') and (indel == 'true') and ( indelcount !=0 ) :
            decine = indelcount*10
            indelcount = decine + int(c)
            indelmax = indelcount
        
        elif c in 'ACGTNacgtn' and indel =='true'   :
            indelcount -= 1
            if indelcount == 0 :
                if insert == 'true' :
                    allelecount['insertions'] += indelmax
                elif insert == 'false' :
                    allelecount['deletions'] += indelmax
                else :
                    print 'ERROR at position',  position
                    
                indelmax = 0
                insert = ''
                indel = 'false'
    
    outfile.write('\t'.join([position, refallele,  coverage]) +'\t' )
    outfile.write('\t'.join(  [ str(allelecount[c] ) for c in allelecountKeys ] ) +'\n' )
    
    cov = int(coverage)
    
    if (cov >=  5 ) and ( allelecount[ref] >= cov/2 )  :          # sequenza identica al riferimento e coverage maggiore (o uguale) a 5 ----> tengo il nucleotide (del riferimento)
        newseq = newseq + ref
        if  cov >=  50 :          # se la coverage e' maggiore (o uguale) a 50  CONTROLLO ANCHE SE HO UNA VARIANTE MINORITARIA (>minallelefreq)
            for c in 'ACGT' :
                if float(allelecount[c])/cov >= minallelefreq  :
                    altc = altc+c
                    altf.append( float(allelecount[c])/cov  )
    
            if altc == ref :
                altc = ''
                altf = []
        
    elif  ( allelecount[ref] < cov/2 ) and (cov >=  20 ) and (cov < 50 ) :       # sequenza diversa dal riferimento e coverage maggiore (o uguale) a 20 e minore di 50 ---->  inserisco una ambiguita'
        myc = ''
        for c in 'ACGT' :
            if allelecount[c] >= allelecount[ref] :
                myc = myc +c 
        myambic= invambidic[myc]
        newseq = newseq + myambic

    elif ( allelecount[ref] < cov/2 ) and (cov >=  50 ) :       # frequenza del riferimento minore del 50% e coverage maggiore (o uguale) a 50 ----> prendo il nucleotide della variante maggiore E CONTROLLO SE HO UNA VARIANTE MINORITARIA (>minallelefreq)
        myc = ref
        for c in 'ACGT' :
            if allelecount[c] > allelecount[myc] :
                myc = c

            if float(allelecount[c])/cov >= minallelefreq :
                altc = altc+c
                altf.append( float(allelecount[c])/cov  )
    
        if altc == myc :
            altc = ''
            altf = []
                
            
        newseq = newseq + myc
        
    else :                                                                  # sequenza diversa dal riferimento e coverage minore di 20 oppure nucleotide qualunque ma coverage minore di 5 ----> metto una "N"
        newseq = newseq + "N"   
    
    positions.append(int(position) )
    
    if altc != '' :
        laltc  = [a for a in altc]
        merged = sorted([ (altf[i],laltc[i] ) for i in range(len(laltc)) ] )
        myaltc = [ el[1]  for el in merged ]
        myalf =  [ el[0]  for el in merged ]
        quasalls.append(myaltc)
        quasfreqs.append(myalf)
        quaspos.append(int(position))

# Adesso scrivo  la sequenza, cercando di evitare che ci siano buchi in mezzo, che sono piu' probabilmente buchi di coverage e non delezioni
# la lunghezza del coronavirus di riferimento NC_045512.2 e' di 29903 nucleotidi :

p=0
mysequence = ['']*29903
for n in range(len(positions)) :
    p=positions[n]
    c=newseq[n]
    mysequence[p-1] = c

lastp=p

init = 0
for i in range(lastp) :
    if (mysequence[i] != '') and  ( init == 0 ) :
        init = 1
    elif (mysequence[i] == '') and  ( init == 1 ) :
        mysequence[i] = 'N'


myLASTsequence = ''.join(mysequence)

outfasta.write(myLASTsequence+ '\n')


#ALLquasispecie = sorted([ (quasfreqs[i], quasalls[i], quaspos[i]  ) for i in range(len(quaspos)) ], reverse=True) 
#
##Q = {}
##
##for q in range( len(ALLquasispecie) ) :
##    myqs = ALLquasispecie[ : ( q+1 ) ]
##    Q[(q+1)] = []
#
#alternatives = []
#alterfreqs = []
#alterpos = []
#
#Q= []
#
#for el in  ALLquasispecie    : # myqs : 
#    alterfreqs.append(el[0])
#    alternatives.append(el[1])
#    alterpos.append(el[2])
#
#protoq = list(itertools.product(*alternatives))
#protofreq = list(itertools.product(*alterfreqs))
#
#for n in range(len(protoq)) :
#    quasispecie = ''.join(protoq[n])
#    quasifreq = numpy.prod(protofreq[n])
#    if quasifreq >= 0.01 :
#        print quasispecie
#        Q.append([quasifreq, quasispecie])   
#        
#
#print "Frequenza minima varianti analizzate:\t" ,  minallelefreq
#print "Numero varianti analizzate:\t",  len(ALLquasispecie)
#print "Numero di quasispecie corrispondenti:\t" ,  len(Q)
#
#n = 0
#for q in Q :
#    n += 1
#    myQsequence = list(mysequence)
#    qfreq = str( round(q[0] , 2) )
#    qseq = q[1]
#    for i in range(len(alterpos)) :
#        myQsequence[ alterpos[i] -1 ] = qseq[i]
#    
#    myQLASTsequence = ''.join(myQsequence)
#    
#    outquas.write(">"+ fastaheader+"_quasisp_"+str(n)+ "_freq_"+ qfreq+ '\n')
#    outquas.write(myQLASTsequence+ '\n')
        

infile.close()
outfile.close()
outfasta.close()
#outquas.close()
