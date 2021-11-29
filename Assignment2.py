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


