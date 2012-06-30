#!/usr/bin/env python
'''Created by Daniel Sindhikara, sindhikara@gmail.com
Program reads an amber input file in proper format and
stores all data in a dictionary which can be looked up
by the FLAG string
'''
import sys
import re

def formatvalues(valuelist,format):
    '''
Returns a string with values in the list 'valuelist'
in format of the string 'format' (e.g. '10I8' or '5e16.8') 
Prepends %FORMAT(....) line
    '''
    valuetype=re.split('\d+',format)[1].lower() # should be e,a,i,...
    numlist=re.split('\D+',format)
    valsperrow=int(numlist[0])
    charspervalue=int(numlist[1])
    myformatstring='%'+"%d" % charspervalue
    formattedstring='%FORMAT('+str(valsperrow)+str(valuetype)+str(charspervalue)
    if len(numlist)>2 : # if there is decimal specification
        numdecimals=int(numlist[2])
        myformatstring += ".%d" % numdecimals
        formattedstring += "."+str(numdecimals)
    
    formattedstring+=")\n"   
    if valuetype == "a" : # Python doesn't understand "a"
        valuetype="s"
    elif valuetype == "I" :
        valuetype="d"
    myformatstring+=valuetype
    remainingonrow=valsperrow
    for value in valuelist:
        remainingonrow-=1
        if remainingonrow<0 :
            formattedstring += "\n"
            remainingonrow=valsperrow-1
        formattedstring+=myformatstring % value

    return(formattedstring)

def readamberfile(filename):
    def convvalfmt(invalues,valuetype):
        #converts string values to appropriate value type
        if valuetype=="e": #exponential
            outvalues = [float(invalue) for invalue in invalues]
        if valuetype=="i": #integer
            outvalues = [int(invalue) for invalue in invalues]
        if valuetype=="a": #character
            outvalues = invalues # no conversion necessary
        return(outvalues)
    amberdic={}
    f = open(filename,"r")
    nextlineisformat=0
    flag = "first"
    values=[]
    for line in f:
        if "FLAG" in line:
            if flag != "first" : # if not the first flag found
                amberdic[flag]=convvalfmt(values,valuetype) # assign previously read value set
                values=[]
            flag=line.split()[1]
            nextlineisformat=1
        elif nextlineisformat==1:
            valuetype=re.split('\d+',line)[1].lower() # should be e,a,i,...
            valspercolumn=int(re.split('\D+',line)[1])
            charspervalue=int(re.split('\D+',line)[2])
            #for now ignore spaces right of decimal
            nextlineisformat=0
        elif "%" not in line:
            valuesinthisline=int(len(line)/charspervalue)
            for i in range(valuesinthisline):
                values.append(line[i*charspervalue:(i+1)*charspervalue]) # add one value

    amberdic[flag]=convvalfmt(values,valuetype) # one last time for last value et
    f.close()
    return(amberdic)

def main():

    minargs=1
    numargs=len(sys.argv)
    if(numargs<minargs) :
    	print "Insufficient arguments, need ",minargs," : ",numargs
        print "format djsamberformat.py <amberfilename>" 
    readamberfile(sys.argv[1])

if __name__ == '__main__' :
    main()
