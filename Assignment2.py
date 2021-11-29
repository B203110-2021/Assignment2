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

