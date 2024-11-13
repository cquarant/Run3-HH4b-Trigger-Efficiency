import os
File_List = open("RootFile_JetMETC_Legacy")
i=1
j=0
fout      = open("RootFile_"+str(j),'w')
fout_name = open("RootFile_Name"+str(j),'w')
for line in File_List:
    fout.write(line)
    fout_name.write(line.split("/")[-4]+"_"+line.split("/")[-1].split("'")[0])
    j+=1
    fout.close()
    fout_name.close()
    fout = open("RootFile_"+str(j),'w')
    fout_name = open("RootFile_"+str(j)+"_Name",'w')
    i+=1
print("queue "+str(j))
File_List.close()