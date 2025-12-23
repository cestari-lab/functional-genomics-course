# Basic Linux commands to get started. 
# Any word after # will not be read as a command line. We use it for notation/information.

# list files
ls	

# Now, try to note the differences between these commands
ls -l 		
ls -lt
ls -ltr
ls -la

# Show the path of your directory
pwd

# mount the drive 
cd /mnt/

#list folders/files
ls -l

# get into your working directory
cd c

# Once in c or d or whatever is your chosen directory, create your working folder
# create folders/directories using mkdir, do not add spaces to your folder names
mkdir myfolder

# open the folder
cd myfolder

# create a text file
nano myfile.txt

#type something "my favourite course is fungenomics!" and close it, saving it.

# see the content using cat
cat myfile.txt

#if you have a lot of content, it is helpful to see the beginning or end of the file
head myfile.txt 	# seeing the beginning few lines of the file
tail myfile.txt	  # seeing the end few lines of the file
less myfile.txt 	# type q to get out of it if you use more or less commands.

# go back to your initial directory
cd ..

# remove (delete) folders and their files
rm -r myfolder

# end your session
exit
