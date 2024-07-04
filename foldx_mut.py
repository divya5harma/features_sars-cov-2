import os
## take mutations and calculate stability and interaction energy

## code is same for paratope as well. minor changes. 

## keep the pdb file and foldX in the same folder. code will automatically create folder for each residues and 2 CSV files with output. Check out for crashes. I suggest to run 4-5 copy of code for different PDBs at the same time.

## list of epitopes
#all epitopes
interface_res = ["A:403:R","A:405:D","A:406:E","A:408:R","A:409:Q","A:415:T","A:416:G","A:417:K","A:420:D","A:421:Y","A:455:L","A:456:F","A:457:R","A:458:K","A:460:N","A:473:Y","A:474:Q","A:475:A","A:476:G","A:477:S","A:486:F","A:487:N","A:489:Y","A:493:Q","A:495:Y","A:502:G","A:505:Y"]

## sometimes foldX will stop in the middle. so you can remove the residues from the list which are done..
## leftovers
interface_res = ["A:420:D","A:421:Y","A:455:L","A:456:F","A:457:R","A:458:K","A:460:N","A:473:Y","A:474:Q","A:475:A","A:476:G","A:477:S","A:486:F","A:487:N","A:489:Y","A:493:Q","A:495:Y","A:502:G","A:505:Y"]

aa = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]

## initialization set chain name and PDB file name.
pdb = "7c01_Repair.pdb"
chain_spike = "A"

# most antibody has H and L chain if this is not the case then change HL in the analyseComplex command.

file4 = open("interaction energy.csv","w")
file5 = open("stability.csv","w")
for i in interface_res:
    file1 = open("individual_list.txt","w")
    temp = i.split(":")
    for j in aa:                                                               ################### change @#################
        file1.write(temp[2]+temp[0]+temp[1]+j+";\n")
    folder = i.replace(":","_")
    file1.close()
    os.system("mkdir "+i.replace(":","_"))
    command1 = ".\\foldx.exe --command=BuildModel --pdb="+pdb+" --mutant-file=individual_list.txt"+" --output-dir "+i.replace(":","_")
    os.system(command1)
    for ct in range(1,21):
        try:
            command2 = ".\\foldx.exe --command=AnalyseComplex --pdb="+pdb.split(".")[0]+"_"+str(ct)+".pdb --analyseComplexChains=HL,"+chain_spike+" --pdb-dir="+i.replace(":","_")+" --output-dir="+i.replace(":","_")
            os.system(command2)
            file2 = open(i.replace(":","_")+"/"+"Summary_"+pdb.split(".")[0]+"_"+str(ct)+"_AC.fxout","r")
            content = file2.readlines()
            
            file2.close()
			
            command3 = ".\\foldx.exe --command=Stability --pdb="+pdb.split(".")[0]+"_"+str(ct)+".pdb --pdb-dir="+i.replace(":","_")+" --output-dir="+i.replace(":","_")
            os.system(command3)
            file3 = open(i.replace(":","_")+"/"+pdb.split(".")[0]+"_"+str(ct)+"_0_ST.fxout","r")
            content2 = file3.readlines()
            file3.close()
            file4.write(i+","+aa[ct-1]+","+content[9].replace("\t",","))
            file5.write(i+","+aa[ct-1]+","+content2[0].replace("\t",","))
        except: ## reducing the chances of failure. it will rerun if it fails.. you can copy it multiple times if you want. Also note that if this second code runs.. it will create duplicate in interaction energy.. make sure nothing repeats in the final output file.
            command2 = ".\\foldx.exe --command=AnalyseComplex --pdb="+pdb.split(".")[0]+"_"+str(ct)+".pdb --analyseComplexChains=HL,"+chain_spike+" --pdb-dir="+i.replace(":","_")+" --output-dir="+i.replace(":","_")
            os.system(command2)
            file2 = open(i.replace(":","_")+"/"+"Summary_"+pdb.split(".")[0]+"_"+str(ct)+"_AC.fxout","r")
            content = file2.readlines()
            file2.close()
			
            command3 = ".\\foldx.exe --command=Stability --pdb="+pdb.split(".")[0]+"_"+str(ct)+".pdb --pdb-dir="+i.replace(":","_")+" --output-dir="+i.replace(":","_")
            os.system(command3)
            file3 = open(i.replace(":","_")+"/"+pdb.split(".")[0]+"_"+str(ct)+"_0_ST.fxout","r")
            content2 = file3.readlines()
            file3.close()
            file4.write(i+","+aa[ct-1]+","+content[9].replace("\t",","))
            file5.write(i+","+aa[ct-1]+","+content2[0].replace("\t",","))

file4.close()
file5.close()

## 2 files are created (interaction and stability).. I suggest to make a master file where you can put all data together. There is some problem in windows that program crashes automatically. so you can rerun for the remaining residues and keep output in one master .xlsx file.
