#    ICRAR - International Centre for Radio Astronomy Research
#    (c) UWA - The University of Western Australia, 2012
#    Copyright by UWA (in the framework of the ICRAR)
#    All rights reserved
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#    MA 02111-1307  USA
#
#
# A.W. [ESO]: 2002-02-28
#             2003-03-18
#             2003-04-16
#             2003-05-21
#             2003-07-04  Fixed bug with COMMENT and HISTORY keywords
#             2003-07-08  CRC32 datasum implemented
#             2003-07-14  Several small changes and fixes
#             2003-08-04  New class structure
#             2003-09-02  Bug fix in fitsHead (calculation of blankCards)
#             2003-12-15  Bug fix in HeadDict in the creation of new header keywords
#             2003-12-15  Bug fix in HeadDict in the creation COMMENT and HISTORY keys
#             2005-06-03  Bug fix in getKeyword and restructure of HeadDict class.
#             2005-07-20  New method FitsHead.parseFitsHead2TupleList and
#                         FitsHead.getKeyType
#             2005-07-21  Added possibility to pass file and cStringIO objects to
#                         FitsHead.__init__
# A.W. [ICRAR]:
#             2019-08-27  Ported to Python3
from printhead import __version__
import sys
import os
import types
import subprocess
import string,re
from glob import glob
from zlib import crc32
from math import ceil

if sys.version_info.major == 3:
    PY_VERSION = 3
else:
    PY_VERSION = 2

if PY_VERSION == 3:
    from io import IOBase

class HeadDict(dict):
    """
    This class defines the data structure used by FitsHead. Essentially the data structure
    consists of nested dictionaries (hash tables) with the following structure:
        {'index':{<index1>:<keyword1>,<index2>:<keyword2>,...},'nodes':{<keyword1>:{'Value':<value1>,'Comment':<comment1>,'Type':<type1>},
                                                                        <keyword2>:{'Value':<value2>,'Comment':<comment2>,'Type':<type2>},...}}

    The 'index' dictionary keeps the sorting information of the keywords, where the key
    is the location of the keyword in the header and the value is the name of the keyword.

    The 'nodes' dictionary contains the keyword name as the key. The value of a single node key
    is again a dictionary. For normal standard keywords the value contains another dictionary,
    (keyval dictionary) which has the three keys 'Value', 'Comment' and 'Type'. For HIERARCH keywords it contains
    the next level in the hierarchy, where the leaf node contains finally a standard keyval
    dictionary as described above for normal keywords.
    """

    def __init__(self, number=0, pos=0):
        """
        """
        self.update({'index':{},'nodes':{}})
        self.POS = pos
        self.NUMBER = number
        self.HEADERSIZE = -1
        self.DATASIZE = (-1,-1)
        self.XmlHead = []


    def setHeaderSize(self, size=-1):
        """
        Set the SIZE variable, which contains the size in bytes of the header.

        INPUT:     size, long
        RETURN:    1 if successful, 0 else
        """
        if PY_VERSION == 3:
            return 1  # Python 3 has a default unlimited integer type
        if type(size) == types.LongType or type(size) == types.IntType:
            self.HEADERSIZE = long(size)
            return 1
        else:
            return 0


    def setPos(self, position=-1):
        """
        Set the POS variable, which contains the position in bytes of the header.

        INPUT:     position, long
        RETURN:    1 if successful, 0 else
        """
        if PY_VERSION == 3:
            return 1  # Python 3 has a default unlimited integer type
        if type(position) == types.LongType or type(position) == types.IntType:
            self.POS = long(position)
            return 1
        else:
            return 0


    def setNumber(self, number=0):
        """
        Set the POS variable, which contains the number of the header counted from 0.

        INPUT:     number, int or long
        RETURN:    1 if successful, 0 else
        """
        if PY_VERSION == 3:
            return 1  # Python 3 has a default unlimited integer type
        if type(number) == types.LongType or type(number) == types.IntType:
            self.NUMBER = int(number)
            return 1
        else:
            return 0


    def setDataSize(self):
        """
        Calculate and set the DATASIZE variable.

        INPUT:     none
        OUTPUT:    int tuple, (datasize, <number of blocks>)
        """
        naxis = int(self.getKeyword('NAXIS')[1])
        siz = 0
        nblocks = 0
        if (naxis > 0):
            siz = 1
            ii  = 1
            while (ii <= naxis):
                kk="NAXIS" + str(ii)
                siz = siz * int(self.getKeyword(kk)[1])
                ii += 1


            siz = siz * abs(int(self.getKeyword('BITPIX')[1]))/8     #calculate data size
            if PY_VERSION == 3:
                nblocks = int(siz/2880.+0.5)
            else:
                nblocks = long(siz/2880.+0.5)

        self.DATASIZE = (siz,nblocks)

        return (siz,nblocks)




    def keyTuple2Dict(self,keyTuple,force=0):
        """
        Method takes a keyTuple on input and returns a dictionary of the HeadDict
        form for that keyword. If the keyword does not exist in self or
        force=1 the dictionary is created from scratch.

        INPUT:     keyTuple, defining a single keyword
        OUTPUT:    keyDictionary
        """

        try:
            (key,value,comment,typ,index) = keyTuple
            hkeys = key.split()

            fullInd = ''
            for hk in hkeys:
                fullInd += "['"+hk+"']"

        except:
            return 0

        existKey = self.getKeyDict(key)
        maxInd = self.getKeyIndex('END')
        newInd = 0

        if list(existKey['index'].keys())[0] == -1:    # keyword not found, create!
            testKey = existKey.copy()
            lenInd = len(self['index'])
            if lenInd > 0:
                maxInd = max(self['index'].keys())

            if maxInd > (lenInd-1) and not lenInd-1 in self['index']:
                newInd = lenInd-1
            elif maxInd > (lenInd-1)  and not len_ind in self['index']:
                newInd = maxInd
            elif maxInd > 0 and self['index'][maxInd] == 'END'  and maxInd == lenInd-1:
                newInd = maxInd
                del(self['index'][maxInd])
                maxInd = maxInd+1
                self['index'].update({maxInd:'END'})
            eval("testKey['nodes']"+fullInd+".update({'Comment':comment,'Value':value,'Type':typ})")
            testKey.update({'index':{newInd:key}})
        elif force == 1:
            testKey = existKey.copy()
            eval("testKey['nodes']"+fullInd+".update({'Comment':comment,'Value':value,'Type':typ})")
        else:
            testKey = existKey


        return testKey


    def getNode(self,key=''):
        """
        Return node of HD dictionary for a certain keyword.

        INPUT:     string, keyword
        OUTPUT:    node of HD dictionary.
        """
        hkeys = key.split()
        fullInd = ''
        node = self['nodes']
        for hk in hkeys:
            node = node[hk]
        return node


    def getElementType(self,key=''):
        """
        Method returns the type of an element ['node'|'leaf'].
        """
        if self.getDescendantElements(key=key) == ('Comment','Type','Value'):
            return 'leaf'
        else:
            return 'node'



    def getDescendantElements(self,key=''):
        """
        Method returns the names of the descendant elements of a <key> as a tuple.

        INPUT:     string, keyword
        OUTPUT:    string tuple, list of descendant nodes
        """
        nodes = self.getNode(key=key).keys()
        return tuple(nodes)



    def getKeyDict(self,key,desc=0,inst=0):
        """
        Method takes a keyword <key> and returns a dictionary of the HeadDict
        form for that key only.
        If desc is different from 0 only the descendant nodes will be returned.
        If inst is different from 0 an HeadDict instance is returned.

        INPUT:     string, keyword
                   int attribute desc, 0 or 1, 0 is default, optional
                   int attribute inst, 0 or 1, 0 is default, optional
        """
        hkeys = key.split()

        fullInd = ''
        curkey = ''
        exists = 0
        keyDict = HeadDict()
        node = self['nodes'].copy()
        for hk in hkeys:
            if hk in node and type(node[hk]) == type({}):
                keys = node[hk].keys()
                exists = 1
                node = node[hk]
                if keys != ('Comment','Value','Type') and hk != hkeys[-1]:
                    if desc == 0:
                        keyDict.getNode(curkey).update({hk:{}})
#                        eval("keyDict['nodes']"+fullInd+".update({hk:{}})")
                elif keys == ('Comment','Value','Type') or hk == hkeys[-1]:
                    if desc == 0:
                        keyDict.getNode(curkey).update({hk:node})
#                        eval("keyDict['nodes']"+fullInd+".update({hk:node})")
                    else:
                        keyDict['nodes'].update(node)

            elif hk in node and type(node[hk]) != type({}):
                node = node[hk]
                if desc == 0:
                    keyDict.getNode(curkey).update({hk:node})
#                    eval("keyDict['nodes']"+fullInd+".update({hk:node})")
                else:
                    keyDict['nodes'].update({hk:node})

            else:
                exists = 0
                if hk == hkeys[-1]:
                    node = {hk:{'Comment':'','Value':'','Type':''}}
                else:
                    node = {hk:{}}
                if desc == 0:
                    keyDict.getNode(curkey).update(node)
#                    eval("keyDict['nodes']"+fullInd+".update(node)")
                else:
                    keyDict['nodes'].update(node)

            fullInd += "['"+hk+"']"
            curkey = (curkey+" "+hk).strip()

        if self.getKeyIndex(key) >=0:
            keyDict['index'].update({self.getKeyIndex(key):key})
        else:
            keyDict['index'].update({-1:key})

        if inst==1:
            return keyDict
        else:
            return keyDict

    def getRegexpKey(self,key):
        """
        Method takes a keyword regular expression string on input and returns
        a list of keywords matching the expression.

        INPUT:     string, keyword regexp
        OUTPUT:    string list, keyword names
        """
        res=[]
        reg=re.compile(key)
        keyList=list(self['index'].values())
        for k in keyList:
            if reg.match(k):
                res.append(k)
        return res



    def getKeyIndex(self,key):
        """
        Return the index of a keyword. Method tests for existence.

        INPUT:    none
        OUTPUT    integer, index of the keyword or -1 if keyword does not exist
        """

        ind = -1
        try:
            list(self['index'].values()).index(key)
            exist = 1
        except:
            exist = 0

        if exist:
            for ind in self['index'].keys():
                if self['index'][ind] == key:
                    return ind

        return ind


    def filter(self,keyexp):
        """
        Remove keywords matching the regular expression <keyexp>

        INPUT:     keyword regular expression
        OUTPUT:    none
        """
        try:
            recomp = re.compile(keyexp)
        except Exception as e:
            return e
        matchlist = map(lambda x,y:(re.match(recomp,x) != None)*(y+1),\
            list(self['index'].values()),\
            list(self['index'].keys()))
        indlist = map(lambda x:x-1,filter(lambda x:x>0, matchlist))
        map(lambda x:self['index'].pop(x), indlist)


    def setKeyIndex(self,newind,key):
        """
        Set the index of a <key> to <newind>.

        INPUT:     int, new index for key
                   string, keyword
        OUTPUT:    1 if succesful, 0 else

        NOTE: This method is pretty destructive, i.e. it really sets the
              index to the keyword given no matter whether this index is
              already occupied or not. It should only be used when creating a
              header from scratch.
        """

        try:
            pos = list(self['index'].values()).index(key)
            self['index'].remove(pos)
            self['index'].update({newind:key})
        except:
            return 0

        return 1


    def deleteKeyIndex(self,key):
        """
        Method takes a keyword <key> and deletes the index in the index
        dictionary of the HeadDict instance. It returns the deleted index.

        INPUT:     string, keyword
        OUTPUT:    int, index of keyword
        """
        ind = self.getKeyIndex(key)
        del(self['index'][ind])
        return ind



    def getKeyword(self,key,check=0):
        """
        Method takes a keyword <key> and returns a tuple of the form
        (<key>,<value>,<comment>,<type>).

        INPUT:  key, string      Name of the keyword to be searched for
        OUTPUT: tuple of strings: (<key>,<value>,<comment>,<type>)

        If the keyword does not exist the strings in the tuple are empty.

        If check is set to 1 and the keyword does not exists, the function
        returns None instead.
        """
        hkeys = key.split()

        fullInd = ''
        for hk in hkeys:
            fullInd += "['"+hk+"']"


        try:
            val = eval("self['nodes']"+fullInd+"['Value']")
        except:
            if check == 1:
                return None
            else:
                val = ''
        try:
            com = eval("self['nodes']"+fullInd+"['Comment']")
        except:
            com = ''
        try:
            if len(str(val)) > 0 and not eval("'Type' in self['nodes']"+fullInd):
                typ = self.getKeyType(key)
            elif eval("'Type' in self['nodes']"+fullInd):
                typ = eval("self['nodes']"+fullInd+"['Type']")
                if typ == '':
                   typ = self.getKeyType(key)
            else:
                typ = ''
        except:
            typ = ''
        return key,val,com,typ,-1



    def updateKeyword(self,keyDict,force=0):
        """
        Method takes a keyDict and updates HeadDict accordingly.

        INPUT:  keyDict, Dictionary of HeadDict structure
                force, optional. If 1 existing keywords are overwritten.
        OUTPUT: 1 for success, 0 otherwise

        """
        keys = list(keyDict['index'].values())
        inds = list(keyDict['index'].keys())

        for ii in range(len(keys)):
            key = keys[ii]
            ind = inds[ii]
            hkeys = key.split()

            fullInd = ''
            curkey = ''

            node = keyDict['nodes']
            oDict = self.getKeyDict(key,inst=1)
            for hk in hkeys:
                if 'Value' in oDict.getNode(key=curkey)[hk]:
                    test = len(oDict.getNode(key=curkey)[hk]['Value'])
                else:
                    test = hk in self.getNode(key=curkey)
                node = node[hk]
                if not test:
                    self.getNode(key=curkey).update({hk:node})
                elif key in ['COMMENT', 'HISTORY', 'ESO-LOG']:
                    if not key in self['nodes']:
                        self['nodes'].update({key:node})
                    self['nodes'][key]['Value'].\
                         append(keyDict['nodes'][key]['Value'][0])
                fullInd += "['"+hk+"']"
                curkey = (curkey+" "+hk).strip()

            if list(self['index'].values()).count(key):
                dind = self.getKeyIndex(key)
                del(self['index'][dind])
            self['index'].update({ind:key})
        return 1


    def getKeyPos(self,key):
        """
        Method takes a keyword <key> and returns the position in the original FITS
        header.

        INPUT:  key, string, Name of the keyword to be searched for
        OUTPUT: int, position of the keyword card in original header

        If the keyword does not exist the output is -1
        """
        values = list(self['index'].values())
        keys = list(self['index'].keys())
        try:
            return keys[values.index(key)]
        except:
            return -1



    def getKeyType(self,key):
        """
        Method updates the keyword dictionary with the derived type
        {<key>:{'Value':<value>,'Comment':<comment>, 'Type':<type>}}
        where <type> is one out of ['B','C','I','L','F','P', 'R' 'T']

        INPUT:     string, keyword
        OUTPUT:    string, derived type or blank string if type could not be derived
        """
        dtRx = re.compile(\
          "(19\d{2}|2\d{3})\-(0\d|1[012])\-([012]\d|3[01])" + \
          "([T ]([01]\d|2[0-3])\:[0-5]\d\:[0-5]\d(\.\d+)?)?\s*$")

        hkeys = key.split()

        fullInd = ''
        for hk in hkeys:
            fullInd += "['"+hk+"']"


        if eval("'Type' in self['nodes']"+fullInd):
            typ = eval("'Type' in self['nodes']"+fullInd)
        else:
            typ = ""
        if eval("'Value' in self['nodes']"+fullInd):

            val = eval("self['nodes']"+fullInd+"['Value']")
            if typ == 'C' or (type(val) == type('') and val.upper() in reserved):
                typ = 'C'
            else:
                try:
                  float(val)
                  value = float(val)
                  if value != 0 and (abs(value) > 1.0e15 or abs(value) < 1e-15):
                      typ = 'R'
                      value = None
                  dotpos = val.find('.')
                  if dotpos < 0 and typ != 'R':
                      try:
                         iv = int(val)
                         if iv < 256 and iv >= 0:
                            typ = 'U'
                         elif iv < 65536:
                            typ = 'S'
                         else:
                            typ = 'I'
                      except:
                         typ = 'L'
                  elif typ != 'R':
                      epos = val.upper().find('E')
                      if epos == -1:
                         epos = len(val)
                      else:
                         ex = int(val[epos+1:])
                      if dotpos >= 0:
                         typ = 'F'
                         if len(val[dotpos+1:epos]) > 15:
                            typ = 'P'
                except:
                    if val == 'T' or val == 'F':
                        typ = 'B'
                    else:
                        typ = 'C'
            if type(val) == type('') and typ == 'C' and dtRx.match(val):
                # check for datetime format
                typ = 'T'


            exec("self['nodes']"+fullInd+".update({'Type':typ})")
            return typ
        else:
            exec("self['nodes']"+fullInd+".update({'Type':''})")

            return ''


    def sortKeys(self):
        """
        Just started a method to do the proper sorting of keywords....
        """

        pk = {'SIMPLE':0,'XTENSION':0,'BITPIX':1,'NAXIS':2,'NAXIS1':3,\
              'NAXIS2':4,'NAXIS3':5,'NAXIS4':6}

        rpk = re.compile('|'.join(pk.keys()))
        keys = list(self['index'].values())
        self['index'] = {}
        keys.sort()
        KeyDict = {}
        kind = -1
        hind = -1
        maxk = -1

        for k in keys:
            if len(k) > 8:
                hind += 1
                KeyDict.update({k:[-2,hind]})
            elif rpk.match(k):
                maxk = max([maxk,pk[k]])
                KeyDict.update({k:[-1,0]})
            elif k != 'END':
                kind += 1
                KeyDict.update({k:[-2,kind]})

        for k in keys:
            if k != 'END':
                if KeyDict[k][0] == -1: KeyDict[k][0] = KeyDict[k][1] + maxk
                if KeyDict[k][0] == -2: KeyDict[k][0] = KeyDict[k][1] + kind + maxk
                self['index'].update({KeyDict[k][0]:k})

        maxk = max(self['index'].keys())
        self['index'].update({maxk+1:'END'})

        return

    def Serialize(self,  dataFl=-1):
        """
        Method creates a list of FITS header cards from the HeadDict dictionary.

        INPUT:     none
        OUTPUT:    string list, FITS header cards

        NOTE: The Dictionary has to be formatted like one of the self
              parts, i.e. self[0]
        """

        specialKeys = ['HIERARCH','COMMENT','HISTORY','ESO-LOG','END']
        self.FHead = []
        FHead = []

        comhist = {'COMMENT':-1,'HISTORY':-1, 'ESO-LOG':-1}

        newind = self['index'].keys()
        newind.sort()

        for ind in newind:

                key = self['index'][ind]
# treat all 'normal' keywords

                if key[0:8] not in specialKeys:
                    fitsLine = key + (8 - len(key))*' '
                    value = str(self['nodes'][key]['Value'])
                    comment = self['nodes'][key]['Comment']
                    typ = self['nodes'][key]['Type']

                    if len(value) > 0:
                        fitsLine += '= '
                        if typ != 'C':
                            fitsLine = fitsLine + (30 - len(fitsLine) - \
                                                   len(value)) * ' '
                        fitsLine = fitsLine + value
                        fitsLine = fitsLine + (39 - len(fitsLine)) * ' '

                    if len(comment) > 0 and len(fitsLine) + 3 < 80:
                        fitsLine = fitsLine + ' / ' + comment

                    if len(fitsLine) > 80: fitsLine = fitsLine[:80]

                    fitsLine = fitsLine + (80 - len(fitsLine)) * ' '
                    if key == 'END':
                        eInd = max(self['index'].keys())+1
                        lInd = len(self)
#                        for ii in range(lInd,eInd):
#                            FHead.append(80*' ')
                        FHead.append((eInd-lInd)*80*' ')
                    FHead.append(fitsLine)

# COMMENT and HISTORY and ESO-LOG keywords

                elif key in ['COMMENT', 'HISTORY', 'ESO-LOG']:
                    comhist[key] += 1
                    if type(self['nodes'][key]['Value']) == type([]):
                        fitsLine = key + self['nodes'][key]['Value'][comhist[key]]
                        fitsLine = fitsLine + (80-len(fitsLine))*' '
                        FHead.append(fitsLine)
                    elif type(self['nodes'][key]['Value']) == type(''):
                        fitsLine = key + self['nodes'][key]['Value']
                        fitsLine = fitsLine + (80-len(fitsLine))*' '
                        FHead.append(fitsLine)



# HIERARCH keywords are reconstructed from the hierarchy in the dict.

                elif key[0:8] == 'HIERARCH':

                    hind = ""
                    hkeys = key.split()

                    fitsLine = key
                    fitsLine = fitsLine + (29 - len(fitsLine)) * ' ' + '= '
                    for hk in hkeys:
                        hind = hind + "['" + hk + "']"

                        kkeys = eval("self['nodes']"+hind+".keys()")

                    if kkeys.count("Value") > 0:
                        value = str(eval("self['nodes']"+hind+"['Value']"))
                        comment = eval("self['nodes']"+hind+"['Comment']")
                        typ = eval("self['nodes']"+hind+"['Type']")
                        if typ != 'C':
                            fitsLine = fitsLine + (43 - len(fitsLine) - \
                                                   len(value)) * ' '
                        fitsLine = fitsLine + value
                        fitsLine = fitsLine + (43 - len(fitsLine)) * ' '
                        if len(fitsLine) + 3 < 80:
                            fitsLine = fitsLine + ' / ' + comment
                        if len(fitsLine) > 80: fitsLine = fitsLine[:80]
                    else:
                        pass

                    fitsLine = fitsLine + (80 - len(fitsLine)) * ' '
                    FHead.append(fitsLine)

                elif key[0:3] == 'END':

                    hlen = len(FHead)
                    blankCards = 36 - ((hlen+1) % 36)
#                    for ii in range(blankCards):
#                        FHead.append(80 * ' ')
                    FHead.append(blankCards*80*' ')
                    FHead.append('END' + 77*' ')

        return FHead




    def XfitsSerialize(self, level=0, indent='   ', pretty=1):
        """
        Method serializes the HD dictionary into a string array. The format is XFits.

        INPUT:     none mandatory
                   int attribute level, >=0 defines the initial indentation level
                                        default 0, optional
                   string attribute indent, whitespace, defines the amount of indentation
                                            per level. default '   ', optional
        OUTPUT:    string list, XML (XFits) formatted header
        """
        if len(self.XmlHead) > 1 and self.XmlHead[0].strip()[:15] == "<HEADER number=":
            return self.XmlHead
        openTags = []
        hflag = 0
        chlist = {'COMMENT':-1,'HISTORY':-1, 'ESO-LOG':-1}
        XmlHead = []
        XmlHead.append(level*indent + '<HEADER number="' + str(self.NUMBER) + '" position="' + \
                           str(self.POS) + '" datasize="' + str(self.DATASIZE[0]) + '">')
        level += 1    # indent all the rest...
        for key in list(self['index'].values()):

# treat all 'normal' keywords
            rkey = key[0:8].strip()
            if rkey not in ['HIERARCH', 'COMMENT', 'HISTORY', 'ESO-LOG']:
                if hflag:
                    for ot in openTags:
                        XmlHead.append(level*indent + '</'+ot+'>')
                        if ot != 'HIERARCH': level = (level - 1) * pretty
                    openTags = []
                    hflag = 0
                if openTags and (openTags[0] in ['COMMENT', 'HISTORY', 'ESO-LOG']):
                    XmlHead.append(level*indent + '</' + openTags[0] + '>')
                openTags = [rkey]
                XmlHead.append(level*indent + '<' + rkey + '>')
                XmlHead.append((level+1)*indent*pretty + '<Value>' + \
                               self.getKeyword(rkey)[1] + '</Value>')
                XmlHead.append((level+1)*indent*pretty + '<Comment>' + \
                               self.getKeyword(rkey)[2] + '</Comment>')
                XmlHead.append(level*indent + '</' + rkey + '>')
                openTags = []

# COMMENT and HISTORY keywords

            elif rkey in ['COMMENT', 'HISTORY', 'ESO-LOG']:
                chlist[rkey] += 1
                if hflag:
                    for ot in openTags:
                        level = (level - 1) * pretty
                        XmlHead.append(level*indent + '</'+ot+'>')
                    openTags = []
                    hflag = 0
                if openTags and openTags[0] != rkey:
                    if len(openTags[0])>0:
                        XmlHead.append(level*indent + '</' + openTags[0] + '>')
                    openTags = [rkey]
                    XmlHead.append(level*indent + '<' + rkey + '>')
                if type(self.getKeyDict(rkey)) == type({}):
                    XmlHead.append(self.getKeyword(rkey)[1])
                elif type(self.getKeyDict(rkey)) == type([]):
                    XmlHead.append(self.getKeyDict(rkey)[chlist[rkey]][1].strip())



# HIERARCH keywords are placed in a real XML hierarchy

            elif rkey == 'HIERARCH':


# COMMENT and HISTORY and ESO-LOG elements are kept open until another element is found
# Close COMMENT or HISTORY here if open
                if openTags and  (openTags[0] in ['COMMENT', 'HISTORY', 'ESO-LOG']):
                    XmlHead.append(level*indent + '</' + openTags[0] + '>')


#hflag controls wheather we are already in a HIERARCH element
#level controls how deep we are in the element and openTags
#keeps all the open tags in reverse order

                if not hflag:
                    level = 2 * pretty
                    XmlHead.append(level * indent + '<HIERARCH>')
                    hflag = 1
                    openTags = ['HIERARCH']

                hind = ""
                hkeys = key.split()
#                    oinds = range(len(openTags))
                oind = 0

#compare current key elements with openTags and close the ones which don't match

                for ot in openTags:
                    if hkeys.count(ot) == 0:
                        XmlHead.append(level*indent + '</'+ot+'>')
                        level -= 1
                        openTags = openTags[1:]
                    else:
                        dum = hkeys.index(ot)
                        oind = max(oind,dum)


                for hk in hkeys[:oind+1]:
                    hind = hind + "['" + hk + "']"


                for hk in hkeys[oind+1:]:
                    level = (level + 1) * pretty
                    hind = hind + "['" + hk + "']"
                    XmlHead.append(level*indent + '<'+hk+'>')
                    openTags = [hk] + openTags

                    kkeys = eval("self['nodes']"+hind+".keys()")
                    if kkeys.count('Value') > 0:
                        dum = kkeys.index("Value")
                        value = eval("self['nodes']"+hind+"['Value']")
                        comment = eval("self['nodes']"+hind+"['Comment']")
                        XmlHead.append((level+1)*indent*pretty + '<Value>' + \
                                       self.getKeyword(key)[1] + '</Value>')
                        XmlHead.append((level+1)*indent*pretty + '<Comment>' + \
                                       self.getKeyword(key)[2]  + '</Comment>')
                        XmlHead.append(level*indent + '</'+hk+'>')
                        level = (level - 1) * pretty
                        openTags = openTags[1:]

# close current <HEADER> element and continue with next
        level -= 1
        XmlHead.append(level*indent + '</HEADER>')


        return XmlHead



    def VotableSerialize(self, level=0, indent='   ', pretty=1):
        """
        Method serializes HeadDict and creates a list of XML strings.
        If <pretty> is 1 (default) then the
        XML file is nicely indented.

        This version is intended to write VOImage output.

        INPUT:     none mandatory
                   int attribute level, >=0 defines the initial indentation level
                                        default 0, optional
                   string attribute indent, whitespace, defines the amount of indentation
                                            per level. default '   ', optional
        OUTPUT:    string list, XML (VOImage) formatted header.
        """

        if len(self.XmlHead) > 1 and self.XmlHead[0].strip()[:13] == "<RESOURCE id=":
            return self.XmlHead
        XmlHead = []
        hflag = 0
        oind = 0
        openTags = ['']

        xstr = '<RESOURCE id="' + str(self.NUMBER) + '"'
        if 'EXTNAME' in self['nodes']:
           xstr = xstr + ' name="' + self['nodes']['EXTNAME']['Value'][1:-1] +'"'

        xstr = xstr + ' type="meta">'
        XmlHead.append(level*indent + xstr)


        level += 1
        XmlHead.append(level*indent + '<INFO name="position" value="' + \
                       str(self.POS) + '"/>')
        XmlHead.append(level*indent + '<INFO name="datasize" value="' + \
                       str(self.DATASIZE[0]) + '"/>')
        for key in self['index'].values():

# treat all 'normal' keywords

#                if key[0:8] != 'HIERARCH':
             if key[0:8] != '--------':     # for test we treat all keywords the same
                 if hflag:
                     for ot in openTags:
                         level = (level-1) * pretty
                         XmlHead.append(level*indent + '</'+ot+'>')
                     openTags = ['']
                     hflag = 0
                 openTags = ['PARAM']

                 (keyword,val,comm,typ,flag) = self.getKeyword(key)

                 if typ == 'I':
                     voTyp = 'int'
                 elif typ == 'U':
                     voTyp = 'unsignedByte'
                 elif typ == 'S':
                     voTyp = 'short'
                 elif typ == 'L':
                     voTyp = 'long'
                 elif typ == 'F':
                     voTyp = 'float'
                 elif typ == 'D':
                     voTyp = 'double'
                 elif typ == 'C':
                     voTyp = 'char'
                 elif typ == 'B':
                     voTyp = 'boolean'
                 else:
                     voTyp = ''

                 if type(val) == type(''):
                    XmlHead.append(level*indent + '<PARAM name="' + \
                                key + '" value="' + val + '" datatype="' +\
                                voTyp + '">')
                    XmlHead.append((level+1)*indent +'<DESCRIPTION>'+\
                                comm + '</DESCRIPTION>')
                    XmlHead.append(level*indent + '</PARAM>')

                 elif type(val) == type([]):
                    for vv in val:
                        XmlHead.append(level*indent + '<PARAM name="' + \
                                    key + '" value="' + vv + '" datatype="' +\
                                    voTyp + '">')
                        XmlHead.append((level+1)*indent +'<DESCRIPTION>'+\
                                    comm + '</DESCRIPTION>')
                        XmlHead.append(level*indent + '</PARAM>')



# HIERARCH keywords are placed in a real XML hierarchy

#                elif key[0:8] == 'HIERARCH':
             elif key[0:8] == '--------':


#hflag controls whether we are already in a HIERARCH element
#level controls how deep we are in the element and openTags
#keeps all the open tags in reverse order

                if not hflag:
                    XmlHead.append(indent + '<PARAM name="' +\
                                   key + '" value="' + self[key]['Value'] + '" datatype="' +\
                                   self[key]['Type'] + '">')
                    hflag = 1
                openTags = ['HIERARCH']
                level = 2 * pretty

                hind = ""
                hkeys = key.split()
                oind = 0

#compare current key elements with openTags and close the ones which don't match

                for ot in openTags:
                    if hkeys.count(ot) == 0:
                        XmlHead.append(level*indent + '</'+ot+'>')
                        level -= 1
                        openTags = openTags[1:]
                    else:
                        dum = hkeys.index(ot)
                        oind = max(oind,dum)


                    for hk in hkeys[:oind+1]:
                        hind = hind + "['" + hk + "']"


                    for hk in hkeys[oind+1:]:
                        level = (level + 1) * pretty
                        hind = hind + "['" + hk + "']"
                        XmlHead.append(level*indent + '<'+hk+'>')
                        openTags = [hk] + openTags

                        kkeys = eval("self"+hind+".keys()")
                        try:
                            dum = kkeys.index("Value")
                            value = eval("self"+hind+"['Value']")
                            comment = eval("self"+hind+"['Comment']")
                            XmlHead.append((level+1)*indent*pretty + '<Value>' + \
                                           value + '</Value>')
                            XmlHead.append((level+1)*indent*pretty + '<Comment>' + \
                                           comment + '</Comment>')
                            XmlHead.append(level*indent + '</'+hk+'>')
                            level = (level - 1) * pretty
                            openTags = openTags[1:]
                        except:
                            pass

# close current <RESOURCE> element and continue with next

        if self.DATASIZE > 0:
            XmlHead.append(level*indent + '<TABLE name="data">')
            level += 1
            XmlHead.append(level*indent + '<FIELD name="image" type="link" '+\
                           'arraysize="[]" datatype="integer">')
            level += 1
            XmlHead.append(level*indent + '<LINK href="cid:' + str(self.NUMBER) + '"/>')
            level -= 1
            XmlHead.append(level*indent + '</FIELD>')
            level -= 1
            XmlHead.append(level*indent + '</TABLE>')

        level -= 1
        XmlHead.append(level*indent + '</RESOURCE>')


# put the XML into the object

        return XmlHead
