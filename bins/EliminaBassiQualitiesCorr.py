
from Bio import SeqIO
import optparse
#import re
#import IUPACfunc




# seleziona solo le reads con qualita' media maggiore di 20 e trasforma il risultato in fastq

parser = optparse.OptionParser(description="QUESTO PROGRAMMA serve ad eliminare le reads con qualita minore di 20 e a salvare in 2 file forward e reverse i formato fasta")
parser.add_option('-i', dest='in_fastq', help='file fastq di input')
parser.add_option('-o', dest='out_fastq', help='nome fastq output: fastq delle sequenze forward e reverse')

parser.add_option('-f', dest='prim_for',  default="NNNNNNNNNN",  help='sequenza primer forward')                                  #
parser.add_option('-r', dest='prim_rev', default="NNNNNNNNNN", help='sequenza primer reverse')                                   #  

parser.add_option('-l', dest='minlength', default="54",  help='lunghezza minima delle sequenze forward e reverse') #<---54 e' la lunghezza dichiarata dell'amplicone piu' corto
parser.add_option('-q', dest='minquality', default="20", help='qualita minima delle sequenze forward e reverse')
parser.add_option('-L', dest='maxlength', default="464",  help='lunghezza massima delle sequenze forward e reverse') #<---464 e' il doppio della lunghezza dichiarata dell'amplicone piu' lungo




args, options = parser.parse_args()

infastq = open(args.in_fastq,'r')
minqual=int(args.minquality)
minlen=int(args.minlength)
maxlen=int(args.maxlength)

primfor=args.prim_for
primrev=args.prim_rev

pforlen=len(primfor)
prevlen=len(primrev)





outfastq=open(args.out_fastq ,  'w')

#outfasta_complement_for =open( args.out_fasta + '_SoloComplement_for.fasta',  'w')

# outfasta_rev= open( args.out_fasta + '_rev.fasta',  'w')

#outfasta_complement_rev = open( args.out_fasta + '_SoloComplement_rev.fasta',  'w')

#outfasta_HighQualities= open(args.out_fasta + '_Scartati_AltaQualita.fasta',  'w')

#outfasta_LowQualities= open( args.out_fasta + '_Scartati_BassaQualita.fasta',  'w')




# TRADUCO GLI AMBIGUITY CODES  #################################################
# ambidic ,  ambinonambidic ,  invambidic ,  complambidic ,  ambikeys  ,  ambiset  = IUPACfunc.iupacdic()

#s = str(primfor)         # ricopio il primer FORWARD com' e' 
#for a in ambikeys:
#    if a in s :
#        s = s.replace(a,  '[' +ambidic[a] + ']' )      
#ambiprimfor=s
#
#s = str(primrev)         # ricopio il primer REVERSE com' e' 
#for a in ambikeys:
#    if a in s :
#        s = s.replace(a,  '[' +ambidic[a] + ']' )      
#ambiprimrev=s

# ##########################################################################


for rec in SeqIO.parse(infastq, "fastq") :
    
    #lunghezza della read:
    mylen=len(rec.seq)
    
    #qualita media della read :
    meanqual = float(sum(rec.letter_annotations["phred_quality"] )  )/mylen
    
    if ( meanqual > minqual ) and (mylen >=minlen ) and  (mylen <= maxlen)  :  #se la qualita' minima e la lunghezza minima e massima sono rispettate
        
        
        SeqIO.write(rec, outfastq, "fastq")    


outfastq.close()

#        
#    #forwards :
#        if  re.match( ambiprimfor, str(rec.seq[:pforlen]) )    : #                (rec.seq[:pforlen] == primfor ) 
#            outfasta_for.write('>'+rec.id+'_FOR_ff\n')
#            outfasta_for.write(str(rec.seq)+'\n')
#
##        elif   re.match( ambiprimfor, str(rec.seq.complement()[:pforlen] )   ) : #                      (rec.seq.complement()[:pforlen] == primfor ) :
##            outfasta_complement_for.write('>'+rec.id+'_FOR_cc\n')
##            outfasta_complement_for.write(str(rec.seq.complement() )+'\n')
#                
##        elif (rec.seq.reverse_complement()[:pforlen] == primfor ) and  (rec.seq.complement()[:prevlen] != primrev ) and (rec.seq[:prevlen] != primrev ) :
##            outfasta_for.write('>'+rec.id+'_FOR_rc\n')
##            outfasta_for.write(str(rec.seq.reverse_complement() )+'\n')
#            
#    #reverse :
#        elif  re.match( ambiprimrev, str(rec.seq[:prevlen]) ) :                                 #           (rec.seq[:prevlen] == primrev ) :
#            outfasta_rev.write('>'+rec.id+'_REV_ff\n')
#            outfasta_rev.write(str(rec.seq)+'\n')              #outfasta_rev.write(str(rec.seq.reverse_complement())+'\n')
#        
##        elif   re.match( ambiprimrev, str(rec.seq.complement()[:prevlen]) )  :        #     (rec.seq.complement()[:prevlen] == primrev )    :
##            outfasta_complement_rev.write('>'+rec.id+'_REV_cc\n')
##            outfasta_complement_rev.write(str(rec.seq.complement() )+'\n')
#                
##        elif (rec.seq.complement()[:prevlen] == primrev ) and  ( rec.seq.reverse_complement()[:pforlen] != primfor ) and  (rec.seq[:pforlen] != primfor ) :
##            outfasta_rev.write('>'+rec.id+'_REV_rc\n')
##            outfasta_rev.write(str(rec.seq.reverse_complement() )+'\n')      
#        
#        #else :
#            #SeqIO.write(rec, outfasta_HighQualities, "fasta")       
#    
#    #else :
#        #SeqIO.write(rec, outfasta_LowQualities, "fasta")       
#            
#
#
#outfasta_for.close()
#
##outfasta_complement_for.close()
#
#outfasta_rev.close()
#
##outfasta_complement_rev.close()
#
##outfasta_HighQualities.close()
##
##outfasta_LowQualities.close()
#
