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


