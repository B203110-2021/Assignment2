#!/usr/local/bin/python3

#import module
import os
import re
import subprocess
import linecache

###########################################################################
###Step1_obtain the relevant protein sequence data and make a preliminary analysis and decision
###########################################################################

print('This is a mini programme that can help you analyze protein sequences. Now, let us get to work!')
#enter working directory name
while True:
 directory=input("please enter your working directory name: ")
#when nothing input, it will warn "Your input is empty!"
 if not directory:
  print('Your input is empty!')
  continue
 else:
#make the working directory
  os.mkdir(directory)
  os.chdir(directory)
  break
 print("All your working result will be put into "+directory)


#enter the protein_family and taxonomic_group that user input
#clear Screen
os.system('clear')
#enter protein family
while True:
 protein_family=input("Enter protine family:\n")
 if not protein_family:
  print('Your input is empty!')
  continue
 else:
  print("Protein family is "+protein_family)
  break

#enter taxonomic group
while True:
 taxonomic_group=input("Enter taxonomic group:\n")
 if not taxonomic_group:
  print('Your input is empty!')
  continue
 else:
  print("Taxonomic group is "+taxonomic_group)
  print("Please wait for a minute because we are getting fasta results!")
  break

#use esearch to search for protein sequences in NCBI database and use efetch to get fasta results
os.system('esearch -db protein -query \"'+protein_family+' in '+taxonomic_group+'\" | efetch -format fasta > protein_seq.fa')
print("All the matched protein sequences in NCBI database are in 'protein_seq.fa'.")

#ask whether user want to check the result
result=input("Do you want to check the result or just continue?\n\tJust type check or continue\n\t")
if result.upper() == "CHECK":
 check=open("protein_seq.fa").read()
 print(check)
else:
 print("So let us continue to check the sequence number!")

#count protein sequence number
count=open("protein_seq.fa").read()
sequence_number=len(list(re.finditer(r'>',count)))

#user's allowable starting sequence set probably shouldn't have more than 1,000 sequences, check if sequence number is more than 1000
#wrong input check
if sequence_number != 0:
#1000 sequence number check
 if sequence_number >= 1000:
  print("The number of sequences is"+" "+str(sequence_number)+" "+"which is more than 1,000 sequences, maybe not suitable to continue analyze!")
 else:
  print("The number of sequences is"+" "+ str(sequence_number) +" "+"which is which is less than 1,000 sequences, so let's do next job!")
else:
#because user's input is wrong, need to quit process and restart 
 print("Your sequence number is 0, your input maybe wrong, please check your protein family name and taxonomic group name carefully and restart the process.")
 quit()

#count species number
#ask the user if they want to continue
q1=input("Do you want to check the number of species or not?\n\tPlease type check or continue\n\t")
if q1.upper() == "CHECK":
 species = open("protein_seq.fa").read()
 re.findall(r"\[.*\]",species)
 species_number=len(set(re.findall(r"\[.*\]",species)))
#show the number of species
 print("The number of species is"+" "+str(species_number))

#ask the user if they want to continue or not continue with the current dataset after checking species number
q2=input("Do you want to continue or not continue with the current dataset after checking species number?\n\tPlease type yes or no\n\t")
if q2.upper() == "YES":
 print("OK!Let us continue our job!")
else:
#because the user do not want to continue with the current dataset after checking species number, so quit process and restart 
 print("OK!We will quit and restart the process.")
 quit()

#do other checks
#count partial sequences
#ask the user if they want to count how many sequences are partial
q3=input("Do you want to check how many sequences are partial or not?\n\tPlease type yes or no\n\t")
if q3.upper() == "YES":
 partial_number=len(list(re.finditer(r'partial',count)))
 print("The number of partial sequences is"+" "+str(partial_number))

#count predicted sequences
#ask the user if they want to count how many sequences are predicted
q4=input("Do you want to check how many sequences are predicted or not?\n\tPlease type yes or no\n\t")
if q4.upper() == "YES":
 predicted_number=len(list(re.finditer(r'PREDICTED',count)))
 print("The number of predicted sequences is"+" "+str(predicted_number))

#count computer automatic annotation protein sequence(start with">XP_**")
#ask the user if they want to count how many sequences are computer automatic annotation
q5=input("Do you want to check how many sequences are automatic annotation or not?\n\tPlease type yes or no\n\t")
if q5.upper() == "YES":
 XP_number=len(list(re.finditer(r'>XP_',count)))
 print("The number of automatic annotation sequences is"+" "+str(XP_number))

###########################################################################
###Step2_determine and plot the level of conservation between the protein sequences
###########################################################################

#ask the user if they want to analyse the level of conservation between the protein sequences
q6=input("Do you want to continue to analyse the level of conservation between the protein sequences or not?\n\tPlease type yes or no\n\t")
if q6.upper() == "YES":
  print("Ok, let us do the analysis!Please wait for the multiple sequence comparisons")
 #Performing multiple sequence comparisons
  subprocess.call('clustalo -i protein_seq.fa -o clustalo.fa --outfmt=clustal -v --force --threads=30',shell=True)
  subprocess.call('clustalo -i protein_seq.fa -o protein.msf --outfmt=clustal -v --force --threads=30',shell=True)

  q7=input("Do you want to see the clustalo result or just continue?\n\tPlease type yes or no\n\t")
  if q7.upper() == "YES":
     os.system('cat clustalo.fa')
     print("These are the clustalo result, please input a window size you want, then please wait for a minute to see the original plotcon.svg on the Firefox.\n\tAfter seeing the plot, we will generate consensus sequences.")
     os.system('plotcon protein.msf -graph svg')
     os.system('firefox plotcon.svg')

#Integration of sequence alignment results to generate consensus sequences
  print("Now let us generate consensus sequences.")
  os.system('cons clustalo.fa clustalo.cons')
  q8=input("Do you want to see the consensus sequences or just continue?\n\tPlease type yes or no\n\t")
  if q8.upper() == "YES":
      con=open("clustalo.cons").read()
      print(con)

#Performing BLASTP
  print("Now let us perform BLASTP")
#Building a protein database
  os.system('makeblastdb -in protein_seq.fa -dbtype prot -parse_seqids -out protein_database')
  print("The protein database has been built, you can see it in your directory.")
  os.system('blastp -db protein_database -query clustalo.cons > blastp_result.out')
  print("BLASTP is done.")
  q9=input("Do you want to see the BLASTP result or just continue?\n\tPlease type yes or no\n\t")
  if q9.upper() == "YES":
       blastp_result=open("blastp_result.out").read()
#only show the summary blastp result
       summary_blastp_result=blastp_result.split(">")
       s_blastp_result=summary_blastp_result[0]
       summary_blast=s_blastp_result[0:-3].split("\n")
       just_line=summary_blast[29:]
       just_line='\n'.join(just_line)
       blast_txt=open('summary_blast.txt','w')
       blast_txt.write(just_line)
       blast_txt.close()
       os.system('cat summary_blast.txt')

#turn the blastp result to list and start filtering
  blast=open('summary_blast.txt').read()
  blast_list=blast.split('\n')
#count the blast result and start filtering
  print("\n\tNow let us count the BLASTP result and start checking and filtering!")
  if len(blast_list)<600:
   print("\n\tThe blast result items number is "+str(len(blast_list))+" ,<600, is profitable,don't need to filter, all the results are in interest_sequence.txt")
#Obtain gene ID 
   ID= []
   for line in blast_list :
#By space to split the sequence
                  ID_line=line.split(' ')
                  ID.append(ID_line[0])
#Line break for gene ID
                  ID_name='\n'.join(ID)
#Writing the gene ID into a file which named interest_sequence.txt
   interest_sequence=open('interest_sequence.txt','w')
   interest_sequence.write(ID_name)
   interest_sequence.close()       
  else:
   print("The blast result items number is "+str(len(blast_list))+" ,>600,the number is a little large,maybe need to filter")
#ask the user if they want to filter the result start with 'XP_' (computer automatic annotation)
   q10=input("Do you want to filter the blast result by delete the result start with 'XP_' (computer automatic annotation)? \n\tPlease type yes or no\n\t")
   if q10.upper() == "YES":
            noXP=[]
            for i in blast_list:
               if re.search(r'^[^(XP_)]',i) :
                 noXP.append(i)
            noXP_content='\n'.join(noXP)
            screen1=open('screen1.txt','w')
            screen1.write(noXP_content)
            screen1.close()
            screen1_r=open('screen1.txt').read()
            print(screen1_r)
            print(len(noXP))
#check the filtered result again
            if len(noXP)<600:
               print("Now The blast result items number <600,is profitable,don't need to filter, all the results are in interest_sequence.txt")
#Obtain gene ID 
               ID= []
               for line in noXP :
#By space to split the sequence
                  ID_line=line.split(' ')
                  ID.append(ID_line[0])
#Line break for gene ID
                  ID_name='\n'.join(ID)
#Writing the gene ID into a file which named interest_sequence.txt
               interest_sequence=open('interest_sequence.txt','w')
               interest_sequence.write(ID_name)
               interest_sequence.close()   
            else:
#If the number of filtered result still >600, just select the first 600 items for follow-up analysis.
               print("Now The blast result items number still >600, we can select the first 600 items for follow-up analysis, the first 600 results are in interest_sequence.txt")
               just_600=[]
               just_600=noXP[0:600]
              #Obtain gene ID 
               ID= []
               for line in just_600 :
#By space to split the sequence
                  ID_line=line.split(' ')
                  ID.append(ID_line[0])
#Line break for gene ID
                  ID_name='\n'.join(ID)
#Writing the gene ID into a file which named interest_sequence.txt
               interest_sequence=open('interest_sequence.txt','w')
               interest_sequence.write(ID_name)
               interest_sequence.close() 

#if the user do not want to filter the result start with 'XP_', ask the user to choose how many sequences he want
   else:
      chosen_seq=[]
      number=input("How many sequences do you want to choose for the analysis? \n\tPlease enter a number which no more than "+str(len(blast_list))+"\n\t")
      chosen_seq=blast_list[0:int(number)-1]
#Obtain chosen_seq ID 
      ID= []
      for line in chosen_seq:
#By space to split the sequence
          ID_line=line.split(' ')
          ID.append(ID_line[0])
#Line break for gene ID
          ID_name='\n'.join(ID)
#Writing the gene ID into a file which named interest_sequence.txt
          interest_sequence=open('interest_sequence.txt','w')
          interest_sequence.write(ID_name)
          interest_sequence.close()   
      print("Now your chosen sequence are in interest_sequence.txt")                  

#Using pullseq to match gene ID with protein sequence
  os.system('/localdisk/data/BPSM/Assignment2/pullseq -i protein_seq.fa -n interest_sequence.txt > interest_sequence.fa')
  q11=input("Do you want to see the interest_sequence.fa or just continue to see the plot of aligned sequences?\n\tPlease type yes or no\n\t")
  if q11.upper() == "YES":
   print(open("interest_sequence.fa").read())

#Using plotcon to get the plot of aligned sequences
  os.system('plotcon -winsize 4 -sformat fasta interest_sequence.fa -graph svg')
  print("The plot which stand the level of conservation between the protein sequences have been save in file plotcon.svg, please wait for a minute and then you can see the plot on Firefox")
  os.system('firefox plotcon.svg')

else:
 print("OK!We will quit and restart the process.")
 quit()


###########################################################################
###Step3_scanning protein sequence(s) of interest with motifs from the PROSITE database
###########################################################################

#ask the user if he want to scan protein sequence of interest with motifs from the PROSITE database
q12=input("Do you want to scan protein sequence of interest with motifs from the PROSITE database or not?\n\tPlease type yes or no\n\t")
if q12.upper() == "YES":
  seperate_fasta=open("interest_sequence.fa").read()
#convert to list
  fasta_list=seperate_fasta.split(">")
#obtain gene ID 
  ID= []
  for x in fasta_list :
    ID_line=x.split(' ')
    ID.append(ID_line[0])

#add the ">" back    
  s='>'
  for a in range(len(fasta_list)):
      fasta_list[a]=s+fasta_list[a]
      new_fasta_list=fasta_list
      

#the list has a superfluous ">" on the first place, delete it
  true_fasta_list=[]
  true_fasta_list=new_fasta_list[1:]
  true_ID=[]
  true_ID=ID[1:]

#convert list to files (the file is named after the gene ID)
  for i in range(len(true_ID)):
        file=open(true_ID[i], 'w')
        file.write(str(true_fasta_list[i]))
        file.close()

#Make a directory whcih called patmatmotifs where user want to put their patmatmotifs results
  os.mkdir("patmatmotifs")
#Using patmatmotifs to scan motifs from the PROSITE database
  Motif_name=[]
  for i in range(len(true_ID)):
          os.system('patmatmotifs -sequence '+true_ID[i]+' -outfile patmatmotifs/'+true_ID[i]+'_motif_result')

  print("All the patmatmotifs results have been put into the directory which named patmatmotifs")     
#generate a summary file that includes all the results, import linecache
#read all files under the path and put them in the list
  root = 'patmatmotifs/'
  file_names = os.listdir(root)
  file_ob_list = []
  for file_name in file_names:
        file_ob = root + file_name
        file_ob_list.append(file_ob)


#for each file, read the contents of the file by line and put it into the same list all_motif
  all_motif = []
  for file_ob in file_ob_list:
        line_num = 1
        length_file = len(open(file_ob, encoding='utf-8').readlines())
        while line_num <= length_file:
            line = linecache.getline(file_ob, line_num)
            line = line.strip()
            all_motif.append(line)
            line_num = line_num + 1

#write the data content to the generated txt file, pay attention to the coding problem
  f = open('./combine_motif_result.txt', 'w+', encoding='utf-8')
  for i, p in enumerate(all_motif):
        f.write(p+'\n')

  f.close()
  print("All the patmatmotifs results have been put into the txt which named combine_motif_result.txt in your directory") 

#ask the user if he want to see all motifs names
q13=input("Do you want to see the result and the names of all motifs or not?\n\tPlease type yes or no\n\t")
if q13.upper() == "YES":
  os.system('cat combine_motif_result.txt')
  Motif_name=[]
  for line in open("combine_motif_result.txt"):
   if 'Motif =' in line:
     Motif_name_line= line.split('Motif = ')[1].split('\n')[0]
     Motif_name.append(Motif_name_line)
print(Motif_name)
print("These are results and the names of all motifs.")


########################################################
#Step4:other appropriate EMBOSS (or other) analysis
########################################################
#module1_sequence display(infoseq)
#ask the user if he want to see the basic information of his inrerest sequence
q14=input("Do you want to see the basic information of your interest sequence?\n\tPlease type yes or no\n\t")
if q14.upper() == "YES":
  os.system('infoseq interest_sequence.fa -outfile interest_sequence.INFO')
  info = open("interest_sequence.INFO").read()
  print(info)

#module2_protein sequence analysis_sequence component statistics(pepstats and compseq)
#ask the user whether he want to do EMBOSS analysis comseq
q15=input("Do you want to do other EMBOSS analysis comseq?\n\tPlease type yes or no\n\t")
if q15.upper() == "YES":
#Make a directory whcih called compseq where user can put their compseq results
 os.mkdir("compseq")
#using compseq to count the occurrence frequency of different strings in sequences according to the specified length
 for x in range(len(ID)):
   os.system('compseq -sequence '+ID[x]+' -out compseq/'+ID[x]+'_compseq_result -word 4')
 print("All the compseq results have been put into the directory which named compseq")

#module2_pepstats
#ask the user whether he want to do EMBOSS analysis pepstats
q16=input("Do you want to do other EMBOSS analysis pepstats?\n\tPlease type yes or no\n\t")
if q16.upper() == "YES":
 os.system('pepstats -sequence interest_sequence.fa -outfile pepstats_result.fa')
 pep=open("pepstats_result.fa").read()
 print(pep)

#module3_secondary structure analysis(garnier and helixturnhelix)
#ask the user whether he want to do EMBOSS analysis garnier
q17=input("Do you want to do other EMBOSS analysis garnier?\n\tPlease type yes or no\n\t")
if q17.upper() == "YES":
 os.system('garnier interest_sequence.fa garnier_result.GARNIER')
 gar=open("garnier_result.GARNIER").read()
 print(gar)

#module3_helixturnhelix
#ask the user whether he want to do EMBOSS analysis like helixturnhelix
q18=input("Do you want to do other EMBOSS analysis helixturnhelix?\n\tPlease type yes or no\n\t")
if q18.upper() == "YES":
 os.system('helixturnhelix -sequence interest_sequence.fa -outfile helixturnhelix_result.fa')
 helix=open("helixturnhelix_result.fa").read()
 print(helix)

#module4_performing BLASTX
#ask the user whether he want to do BLASTX
q19=input("Do you want to do BLASTX?\n\tPlease type yes or no\n\t")
if q19.upper() == "YES":
 query_fasta=input("Enter query_fasta:\n")
 file=open('query_fasta.fasta','w')
 file.write(query_fasta)
 os.system('blastx -db protein_database -query query_fasta.fasta -outfmt 7 > blastx_result.out')
 bx=open("blastx_result.out").read()
 print(bx)

#The end
print("All analysis result have been put into your directory, I hope this programme is useful. Good bye~")



