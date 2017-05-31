# -*- coding: utf-8 -*-
'''
Created on Mon Oct 19 17:30:27 2015
Provides some help for nipype workflows 
(to be used within connect statements mostly)
@author: Pierre-Yves Herve
'''

def getElementFromList(inlist,idx,slc=None):
    '''
    For selecting a particular element or slice from a list 
    within a nipype connect statement.
    If the slice is longer than the list, this returns the list
    '''
    if not slc:
        outlist=inlist[idx]
    else:
        outlist=inlist[idx:slc]
    return outlist

def getElementFromListofList(inlist,idx,idx2):
    '''
    For selecting a particular element from a list of lists 
    within a nipype connect statement.
    '''
    return inlist[idx][idx2]
    
def reverseOrder(inlist):
    '''
    Reverses the list
    '''
    inlist.reverse()
    return inlist

def toList(x):
    '''
    We can't put a built-in inside a connect statement
    Converts to list
    '''
    return [x]

def toListOfLists(mylist):
    '''
    Transforms each element of a list in a list
    Useful for MRConvert interface in a MapNode
    '''
    return [[x] for x in mylist]

def appendToList(mylist,element):
    '''
    Appends element to list
    '''
    if type(mylist) != list:
        # sanity check
        mylist=[mylist]
    mylist.append(element)
    return mylist

def appendToListFlatten(mylist,element):
    '''
    Appends element to list
    '''
    flat =[]
    if type(mylist) != list:
        # sanity check
        mylist=[mylist]
    for elt in element:
        flat += elt
    mylist.append(flat)
    return mylist

# this one for a 3dTShift oddity
def prependArobase(inputstr):
    '''
    Prepends an arobase to string (AFNI related)
    '''
    inputstr="@"+inputstr
    return inputstr
# this one is really really dumb, but sometimes its the only way I find to get merge nodes working...
def getList(in_list):
    '''
    Tautology (for merge nodes)
    '''
    out_list = in_list
    return out_list
  
def concatenateLists(in_lists):
    '''
    Concatenates list inputs
    '''
    return [a + a for a in in_lists]

def unique(in_list):
    '''
    remove duplicates from the list
    '''
    return list(set(in_list))  

def uniqueUnlist(in_list):
    '''
    remove duplicates from the list, 
    and return a single value if all values are the same.
    '''
    unq = list(set(in_list))
    if len(unq)==1:
        return unq[0]
    else:
        return unq

def createSimpleList(*args):
    '''
    makes a list from supplied arguments
    '''
    return(args)

def createList(**kwargs):
    '''
    Puts the keyworded arguments into a list, 
    or concatenate lists if lists are supplied.
    Use when order is not an issue!
    '''
    out_list=[]
    for key, value in kwargs.iteritems():
        if not type(value)==list:
            out_list.append(value)
        else:
            out_list=out_list+value
    return(out_list)

def createList2(item1,item2):
    '''
    Creates a list from the 2 arguments
    This preserves the order in which the arguments are supplied
    '''
    return([item1,item2])

def createList3(item1,item2,item3):
    '''
    Creates a list from the 3 arguments
    This preserves the order in which the arguments are supplied
    '''
    return([item1,item2,item3])


def createList4(item1,item2,item3,item4):
    '''
    Creates a list from the 4 arguments
    This preserves the order in which the arguments are supplied
    '''
    return([item1,item2,item3,item4])

def createList5(item1,item2,item3,item4,item5):
    '''
    Creates a list from the 5 arguments
    This preserves the order in which the arguments are supplied
    '''
    return([item1,item2,item3,item4,item5])

def createList6(item1,item2,item3,item4,item5,item6):
    '''
    Creates a list from the 6 arguments
    This preserves the order in which the arguments are supplied
    '''
    return([item1,item2,item3,item4,item5,item6])


def revString(string):
    '''
    reverses the string
    '''
    string=string[::-1]
    return(string)
    

def build_deformation_series(series1,series2):
    '''
    ANTS deformation series need to be supplied in
    a particular order.
    '''
    series1.reverse()
    series2.reverse()
    deformation_series=series2+series1
    return deformation_series

# to enable meta interfaces
def writeIdentity(**kwargs):
    '''
    Tautology with a dictionary
    Has its uses for meta-interfaces
    '''
    return kwargs# other dictionary functions

def evalIdentity(**kwargs):
    '''
    Evals a python file, outputs a new dict with the 
    value in the file associated with each key
    '''
    return { key:eval(value) for key,value in kwargs.iteritems() }

def getFromDict(in_dict, key):
    '''
    Filters dict input for a particular key
    '''
    return in_dict['key']
    
def evalFile(in_file):
    '''
    evaluates the inout and returns the input
    '''
    with open(in_file,'r') as f:
        out=f.read()
    try:
        out_value=eval(out)
    except:
        out_value=str(out)   
    return out_value

