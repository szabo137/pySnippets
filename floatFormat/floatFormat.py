"""
format a float
"""


def formatF(floatString):
    tempPre,tempPost = floatString.split('.')


    if int(tempPost)>100:
        if tempPost=="9"*len(tempPost):
            tempPost="-10^{%s}"%(len(tempPost))
            tempPre=str(int(tempPre)+1)
    else:
        tempPost = "." + tempPost

    return tempPre + tempPost


if __name__=='__main__':
    testFloat = str(0.9993)
    print formatF(testFloat)
