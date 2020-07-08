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
from printhead.classes.HeadDict import HeadDict

if sys.version_info.major == 3:
    PY_VERSION = 3
else:
    PY_VERSION = 2

if PY_VERSION == 3:
    from io import IOBase

class FitsHead:
    """
    Class parses headers of FITS files and creates a memory data structure
    or creates header files. It supports compressed files
    (.gz and .Z)
    for more details just call the usage function or run the
    script without parameters.
    """
    def __init__(self,file,skey='END',struct=0,show=0,check=0, verbose=0, mode=1):
        """
        """
        self.verbose = int(verbose)
        self.nbytes = 0             # number of bytes read so far
        self.POS = []                # position of headers
        self.SIZE = []
        self.datasum = []            # datasum of headers if check!=0
        self.show = int(show)        # print the header if show!=0
        self.struct = int(struct)    # examine the structure of the file
        self.check = int(check)      # calculate datasums
        self.Extension = []          # list of HeadDict instances
        self.Mode = mode             # if 0 it is assumed that the input does
                                     # not contain data (.hdr file)
        self.KKeys = ['SIMPLE','EXTEND','NAXIS[0-9]{0,2}','BITPIX','XTENSION', 'END',]
        if skey != 'END': self.KKeys.append(skey)
        if type(file) == type(''):
            (self.fd, self.size) = self.openFile(file)
            if self.size == -1:
                errMsg = "*** File %s does not exists ****" % file
                raise Exception(errMsg)
        elif PY_VERSION == 2 and type(file) == types.FileType:  # if fd != 0 we assume that a file object is passed
            (self.fd, self.size) = (file, file.tell())
            self.ID = self.fd.name
            self.name = self.fd.name
        elif PY_VERSION == 3 and isinstance(file, IOBase):  # if fd != 0 we assume that a file object is passed
            (self.fd, self.size) = (file, file.tell())
            self.ID = self.fd.name
            self.name = self.fd.name
        elif type(file).__name__ == 'StringI': # this is a cStringIO object
            (self.fd, self.size) = (file, file.tell())
            self.ID = ""
            self.name = ""
        else:
            errMsg = "Invalid type passed to file parameter during __init__"
            raise Exception(errMsg)
        self.HEAD = []               # list of list(s) of header cards
        self.analyzeStruct()

    def analyzeStruct(self):
        """
        Method does minimal parsing of the headers in order to derive the structure
        of the FITS file. It fills the string array self.STRUCT and the string array
        self.HEAD, which contains the plain header cards. The mandatory keywords
        are parsed into the HD dictionaries for each extension.
        """
        self.STRUCT = []
        HH = self.dumpHead()
        hcount = 1
        headfl = 1
        if self.struct > 0:
            while len(HH) > 0 :
                self.HEAD.append(HH)
                if self.Mode:
                    self.skipData(header=-1)
                naxis = int(self.Extension[-1].getKeyword('NAXIS')[1])
                if headfl == 1:
                    stmp = "# HDR  NAXIS  "
                    for na in range(1,naxis+1):
                        stmp += "NAXIS%d  " % na
                    stmp += '        POS         DATASUM'
                    self.STRUCT.append(stmp)
                    self.STRUCT.append(70*'-')
                    headfl = 0
                if self.check:
                    datasum = self.datasum[-1]
                else:
                    datasum = -1
                stmp = "%3d  %3d    " % (len(self.HEAD), naxis)
                for na in range(1,naxis+1):
                    lna = int(self.Extension[-1].getKeyword('NAXIS'+str(na))[1])
                    stmp += "%6d   " % lna
                if naxis > 0:
                    stmp += "%10d    %12d" % (self.POS[-1][0],datasum)
                self.STRUCT.append(stmp)
                if self.show == len(self.HEAD)-1 and self.show != 99:
                    break
                else:
                    HH = self.dumpHead()
                    hcount += 1
        else:
            self.HEAD = [HH]



    def dumpHead(self):
        """
        Read all header blocks starting at current position.

        Output: HEAD: array of header strings
        """

        _BLOCKSIZE_ = 2880
        skey = self.KKeys[-1]
        rkkeys = self.KKeys[0]
        for kkey in self.KKeys[1:]:
            rkkeys = rkkeys + '|' + kkey
        rq = re.compile(rkkeys)

        index = 0
        number = len(self.Extension)
        endfl = 0
        skfl = 0
        keys=[]
        block = self.fd.read(_BLOCKSIZE_)
        block = block.decode("latin-1")
        self.nbytes = self.nbytes + _BLOCKSIZE_
        if len(block) > 0 and not block[0:8] == 'XTENSION' and not block[0:6] == 'SIMPLE':
            return ''
        if block:
            self.POS.append([self.nbytes - _BLOCKSIZE_,0])
            HD = HeadDict(number=number, pos = self.nbytes - _BLOCKSIZE_)
        HEAD = block
        sline = ''
        while block:
            kkeys=[]
            for ind in range(0,_BLOCKSIZE_,80):
                if block[ind] != ' ':
                    pkey = block[ind:ind+8].strip()
                    if pkey == 'END':
                        endfl = 1
                        key = pkey
                    elif pkey == 'HIERARCH':
                        eqind = block[ind:ind+80].find('=')
                        key = block[ind:ind+eqind].strip()
                    else:
                        key = pkey
                    kkeys.append(key)
                    if rq.match(key):
                        LineTuple = self.parseFitsCard(block[ind:ind+80])
                        sline = block[ind:ind+80].strip()
                        if skey != 'END' and LineTuple[0] == skey:
                            HEAD = sline
                            skfl = 1
                        LineDict = HD.keyTuple2Dict(LineTuple)
                        LineDict['index']={index:LineTuple[0]}
                        HD.updateKeyword(LineDict)
                index += 1

            keys.append(kkeys)
            if endfl == 1:
#               stat=self.fd.close()
                break
            block=self.fd.read(_BLOCKSIZE_)
            block = block.decode("latin-1")
            self.nbytes = self.nbytes + _BLOCKSIZE_
            if skfl == 0: HEAD = HEAD + block


        if block or index > 0:
            HD.setHeaderSize(self.nbytes - self.POS[-1][0])
            HD.setDataSize()
            self.Extension.append(HD)
            self.POS[-1][1] = self.nbytes

        return HEAD


    def skipData(self,header=-1):
        """
        skipData method for multiple extension files. Contains also the calculation of the
        data checksum. If the file object is not created from a ordinary file, like a socket or
        a pipe then the method does not skip but rather read through the data.
        """
        (siz,nblocks) = self.Extension[header].DATASIZE
        siz = int(siz)
        rr = siz % 2880
        checksum = -1
        if (siz > 0):
            if dir(self.fd).count('name') != 0 and (not self.check) and \
                self.fd.name[1:-1] != 'fdopen':    #this fd.name means pipe, i.e. no seek
                if siz != 0: self.fd.seek(siz,1)     #skip over data
                if rr  != 0: self.fd.seek(2880-rr,1) #and rest of card
            else:
                datasiz = siz
                if rr!=0: datasiz = datasiz + (2880-rr)
                data = self.fd.read(datasiz)
                checksum = -1
                checksum = crc32(data)

            self.nbytes = self.nbytes + siz
            if rr != 0: self.nbytes = self.nbytes+(2880-rr)
        else:
            datasiz = 0
            checksum = -1

        self.datasum.append(checksum)
        self.SIZE.append(siz)
        return 0



    def getData(self,header=0,ofile='',blfl=1):
        """
        Method reads and optionally writes the data part at the current position
        in a FITS file.
        If blfl is 0 the actual data size as given in the header is read, else
        the number of complete FITS blocks are read.
        """
        if self.size>0:   # positioning does not work for streams
            self.fd.seek(self.POS[header][0]+self.POS[header][1],0)
        else:
            header = -1   # force header to be last one
        (siz,nblocks) = self.Extension[header].DATASIZE
        wfl = 0
        if len(ofile) > 0:
            try:
                of = open(ofile,'w')
                wfl = 1
            except:
                print("Problem opening output file:",ofile)
                return
        if wfl:
            for ii in range(nblocks):
                block = self.fd.read(2880)
                of.write(block)
            del(block)
            return -1
        else:
            if blfl == 0:
                rsiz = siz
            else:
                rsiz = nblocks*2880
            data = self.fd.read(rsiz)
            return data


    def openFile(self,file):
        """
        Opens the file or a pipe if the file is compressed and returns
        a file-descriptor and the size of the file.
        """
        flist = glob(file)        #try to find the file
        if len(flist) == 0:            # don't open new one if it does not exist
            return (-1,-1)
        else:
            base = os.path.basename(file)
            ID, ext = os.path.splitext(base)
            if ext == '.Z' or ext == '.gz':
                (fd,STDIN,STDOUT) = subprocess.Popen(['/usr/bin/gunzip','-qc',file],0)
                size = -2   # size is not available in a pipe, but this is not a problem
                self.name = file
                self.ID, ext = os.path.splitext(ID)
            else:
                if PY_VERSION == 3:
                    fd=open(file, 'rb')
                else:
                    fd=open(file)
                fd.seek(0,2)
                size = fd.tell()
                fd.seek(0,0)
                self.name = fd.name
                self.ID = base

        return (fd,size)



    def parseFitsHead(self):

        """
        Method parses self.HEAD into a HeadDict dictionary.
        """
        exts = []
        for ii in range(len(self.HEAD)):
            HD = HeadDict()
            for ind in range(0,len(self.HEAD[ii]),80):
                h = self.HEAD[ii][ind:ind+80]
                LineTuple = self.parseFitsCard(h)
                key = LineTuple[0]
                if key in ['COMMENT', 'HISTORY', 'ESO-LOG']:
                    LineDict = {'index':-1,'cards':{key:{'Value':LineTuple[1],\
                                'Comment':LineTuple[2],'Type':''}}}
                else:
                    LineDict = HD.keyTuple2Dict(LineTuple)
                if len(key) > 0:
                    LineDict.update({'index':{ind/80:LineTuple[0]}})
                    HD.updateKeyword(LineDict)

            HD.setNumber(ii)
            HD.setPos(self.Extension[ii].POS)
            HD.setDataSize()
            exts.append(HD)
        self.Extension = exts
        return


    def parseFitsHead2TupleList(self, forceString = 1):

        """
        Method parses self.HEAD into a list containing
        tuples of the form (fileId, ext_ind, key_ind, key, value, comment, type).
        If abs(forceString) == 2 then the DBCM format is produced which contains
        in addition for the numeric types a kw_value_numeric and for certain
        keywords a kw_value_datetime.

        INPUT: forceString    int, optional parameter, if > 0 all entries are converted to strings
                              if < 0 the HISTORY and COMMENT and ESO-LOG keys are not converted.
                              if abs(forceString) == 1: standard format
                              if abs(forceString) == 2: DBCM format
        RETURN: list, list of line tuples.
        """
        if forceString == 0: forceString = 1
        # dateTimeKeys = ['DATE', 'DATE-OBS', 'HIERARCH ESO OBS START', 'HIERARCH ESO TPL START', \
        #                'HIERARCH ESO TEL DATE', 'HIERARCH ESO INS DATE']
        tupleList = []
        for ii in range(len(self.HEAD)):
            tupleList.append([])
            for ind in range(0,len(self.HEAD[ii]),80):
                h = self.HEAD[ii][ind:ind+80]
                LineTuple = self.parseFitsCard(h, index=ind/80)
                key = LineTuple[0]
                LineList = []
                if len(key) > 0:
                    if key in ['COMMENT', 'HISTORY', 'ESO-LOG']:
                        LineList = [self.ID, str(ii), str(LineTuple[4]), LineTuple[0], LineTuple[1][0],\
                             '',LineTuple[3]]
                    else:
                        LineList = [self.ID, str(ii), str(LineTuple[4]), LineTuple[0], LineTuple[1], LineTuple[2], \
                        LineTuple[3]]
                    if abs(forceString) == 2:
                        if LineTuple[3] not in ['C','B','R','T']:
                            kw_value_numeric = LineTuple[1]
                        else:
                            kw_value_numeric = ''
                        if LineTuple[3] == 'T':
                            kw_value_datetime = LineTuple[1][:23].replace('T', ' ')
                        else:
                            kw_value_datetime = ''
                        dotPos = self.ID[2:].find('.') + 2    # make sure that a '.' in the first two characters is ignored
                        if dotPos > 10: dotPos = 10           # and limit the prefix to the first 10 characters
                        LineList = [self.ID[:dotPos]] + LineList + [kw_value_numeric, kw_value_datetime]
                    if forceString > 0 or (key not in ['COMMENT', 'HISTORY', 'ESO-LOG']):
                            tupleList[ii].append(tuple(LineList))
        return tupleList


    def parseFitsCard(self,line, index=-1):
        """
        Method to parse a single FITS header card.

        INPUT: string(80), One line of a FITS header
        RETURN: tuple, (key, value, comment, type, index)
        """

        key = ''
        value = ''
        comment = ''
        typ = ''
        sexpr = re.compile('^COMMENT|HISTORY|END|ESO-LOG')  # these are the special keywords
        qexpr = re.compile("'(''|[^'])*'")  # this allows to catch crazy keyword values like "'o''Neill'"


        if line[0] != ' ' and not sexpr.match(line):
            (key,rest) = line.split('=',1)
            key = key.strip()
            rest = rest.strip()
            if rest[0] == "'":
                try:
                    m = qexpr.match(rest)
                    value = m.group()[1:-1].strip()
                    vind = m.end()
                    typ = 'C'
                except Exception:
                    errMsg = "Could not match expression '(''|[^'])*' against %s\n" % (rest)
                    errMsg += "FITS value of type string not properly quoted! Suspect card is:\n"
                    errMsg += line
                    raise Exception(errMsg)
                if rest[vind:].find("/") > -1:
                    comment = rest[vind:].split("/",1)[1]
                else:
                    comment = ''
            else:
                if rest.find("/") > -1:
                    (value,comment) = rest.split("/",1)
                else:
                    comment = ''
                    value = rest

            value = value.strip()
            comment = comment.strip()
        elif sexpr.match(line):
            key = sexpr.match(line).group()
            rest = sexpr.split(line)
            comment = ''
            if key in ['COMMENT', 'HISTORY', 'ESO-LOG']:
                value = [rest[1].strip()]
            else:
                value = ''

        return self.getKeyType((key,value,comment,typ,index))




    def xmlHead(self, format='vo', outfile = '', pretty = 1, head=0):
        """
        Method takes Extension and creates a list of XML strings or writes
        the XML strings to <outfile>. If <pretty> is 1 (default) then the
        XML file is nicely indented.
        """

        if head == 99:
            heads = self.Extension
        else:
            heads = [self.Extension[head]]
        XmlHead = []

        level = 0
        indent = '   '
        XmlHead.append('<?xml version="1.0" encoding="ISO-8859-1"?>')

        if format == 'vo':

            XmlHead.append('<VOTABLE version="1.1">')
            level = 1 * pretty
            XmlHead.append(level*indent + '<INFO name="Creator" value=' + \
                           '"ESO printhead tool"/>')
            XmlHead.append(level*indent + '<INFO name="Version" value="' +\
                           __version__ + '"/>')
            XmlHead.append(level*indent + '<INFO name="Compatibility" value=' +\
                           '"FITS"/>')
            XmlHead.append(level*indent + '<DESCRIPTION>')
            level +=1

            XmlHead.append(level*indent + 'VOTable file created from FITS file')
            XmlHead.append(level*indent + self.name)

            level -=1
            XmlHead.append(level*indent + '</DESCRIPTION>')
        elif format == 'xfits':
            level = 1 * pretty
            XmlHead.append('<?xml-stylesheet type="text/xml" href="XMLmenu.xsl"?>')
            XmlHead.append('<XFits>')
        else:
            XmlHead.append("<ERROR>Invalid format specified. Should be vo or xfits</ERROR>")
            return XmlHead
        for HD in heads:
            if format == 'vo':
                XmlHead.append(HD.VotableSerialize(level=level, indent=indent))
            else:
                XmlHead.append(HD.XfitsSerialize(level=level, indent=indent))

# close the root element

        if format == 'vo':
            XmlHead.append('\n</VOTABLE>')
        else:
            XmlHead.append('\n</XFits>')


# put the XML into the object

        return XmlHead






    def walkDict(self,level,dict):
        """
        Helper method to recursivly walk-through a dictionary and serialize it.
        """

        for k in dict.keys():
                if type(dict[k]) == type({}):
                    level += 1
                    self.XHead.append((level*"   ") + "<" + k + ">")
                    self.walkDict(level,dict[k])
                    self.XHead.append((level*"   ") + "</" + k + ">")
                else:
                    level += 1
                    if len(dict[k]) > 0:
                        if dict[k][0] == "'" or k == 'Comment':
                            xtype = ' type="string"'
                            self.XHead.append((level*"   ") + "<" + k + xtype + ">" + \
                                      dict[k][1:-1] + "</" + k + ">")
                        else:
                            xtype = ' type="numeric"'
                            self.XHead.append((level*"   ") + "<" + k + xtype + ">" + \
                                              dict[k] + "</" + k + ">")
                    else:
                            self.XHead.append((level*"   ") + "<" + k + "/>")



                level -= 1
        return


    def xmlHead_WalkDict(self,outfile = ''):
        """
        Method takes HeadDict and creates a list of XML strings or writes
        the XML strings to <outfile>.
        This version does not keep the internal order of the HIERARCH keywords, but
        is a lot nicer code than xmlHead and uses the walkDict method recursively.
        """

        self.XHead = []

        self.XHead.append('<?xml version="1.0" encoding="ISO-8859-1"?>')
        self.XHead.append('<XFits>')
        hflag = 0
        for HD in self.Extension:
            for ii in range(len(self)):
                key = self[ii][0:8]
                if key not in ['HIERARCH','COMMENT','HISTORY', 'ESO-LOG']:
                    self.XHead.append('   <' + key + '>')
#                    hkeys = self.walkDict(1,HD[key])
                    self.XHead.append('   </' + key + '>')

                elif key == 'HIERARCH' and hflag == 0:
                    hflag = 1
                    self.XHead.append('   <HIERARCH>')
#                    hkeys = self.walkDict(1,HD['HIERARCH'])
                    self.XHead.append('   </HIERARCH>')

                elif key == 'COMMENT':
                    self.XHead.append('   <COMMENT>')
                    self.XHead = self.XHead + HD['COMMENT']
                    self.XHead.append('   </COMMENT>')

                elif key == 'HISTORY':
                    self.XHead.append('   <HISTORY>')
                    self.XHead = self.XHead + HD['HISTORY']
                    self.XHead.append('   </HISTORY>')
                elif key == 'ESO-LOG':
                    self.XHead.append('   <ESO-LOG>')
                    self.XHead = self.XHead + HD['ESO-LOG']
                    self.XHead.append('   </ESO-LOG>')


        self.XHead.append('</XFits>')

        if len(outfile) > 0:
            try:
                o = open(outfile,'w')
            except:
                print("ERROR: Unable to open ",outfile)
                return 1

            for xml in self.XHead:
                if type(xml) == type(''):
                    o.write(xml + "\n")
                elif type(xml) == type([]):
                    for x in xml:
                        o.write(x + "\n")

            o.close()
        return 0


    def getAllFHeads(self):
        """
        Just for convenience. Not very useful
        """
        self.FHead = []
        for HD in self.Extension:
            self.FHead.append(HD.Serialize())

        return



    def getFitsCard(self,key,header=0):
        """
        Method takes a keyword <key> and returns the original FITS
        card.
        The value of the header parameter can be used to select a
        specific header (counted from 0).

        INPUT:  key, string      Name of the keyword to be searched for
                header,int       number of header to be searched (def. 0)

        OUTPUT: string

        If the keyword does not exist the output string is empty.
        """
        ind = self.Extension[header].getKeyPos(key)
        if ind > -1:
            return self.HEAD[header][ind*80:ind*80+80]
        else:
            return ''


    def getKeyType(self,lineTuple):
        """
        Method tries to guess the type of a keyword value,
        where <type> is one out of ['B','C','U', 'S', 'I','L','F','P', 'R','T']
        and updates the lineTuple on input.

        types:
            'B':    boolean
            'C':    character
            'T':    tiny int unsigned (0 >= value < 256)
            'S':    short (-65536 < value < +65536)
            'I':    integer (32 bit int)
            'L':    long (unlimited)
            'F':    float (if len(value) <= 10)
            'D':    double (if len(value) > 10)
            'T':    datetime string (ISO)

        INPUT: tuple, lineTuple
        RETURN: lineTuple
        """
        # regexp for dateTime type
        dtRx = re.compile(\
          "(19\d{2}|2\d{3})\-(0\d|1[012])\-([012]\d|3[01])" + \
          "([T ]([01]\d|2[0-3])\:[0-5]\d\:[0-5]\d(\.\d+)?)?\s*$")

        # deal with reserved words, which would lead to the wrong type...
        reserved = ['INFINITY', 'INF', 'NAN']
        val = lineTuple[1]
        typ = lineTuple[3]

        # check if type is defined already or if one of the reserved words is used.
        if typ == 'C' or (type(val) == type('') and val.upper() in reserved):
            typ = 'C'
            value = val
        else:
            try:
              float(val)
              value = float(val)
              if value != 0 and (abs(value) > 1.0e15 or value < 1e-15):
                  typ = 'R'
                  value = None
              dotpos = val.find('.')
              if dotpos < 0 and typ != 'R':
                  try:
                     iv = int(val)
                     if iv < 256 and iv >= 0:
                        typ = 'U'
                     elif iv > iv < 65536:
                        typ = 'S'
                     elif iv :
                        typ = 'I'
                     value = iv
                  except:
                     typ = 'L'
                     value = int(val)
              elif typ != 'R':
                  epos = val.upper().find('E')
                  if epos == -1:
                     epos = len(val)
                  else:
                     ex = int(val[epos+1:])
                  if dotpos >= 0:
                     typ = 'F'
                     if len(val[dotpos+1:epos]) > 15:
                        typ = 'D'
            except:
                if val == 'T' or val == 'F':
                    value = True if val=='T' else False
                    typ = 'B'
                else:
                    typ = 'C'
                    value = val
        if type(val) == type('') and typ == 'C' and dtRx.match(val):
            # check for datetime format
            typ = 'T'

        return (lineTuple[0], value, lineTuple[2], typ, lineTuple[4])




    def setVerbose(self, value):
        """
        Set verbosity of output.
        """
        self.verbose = value

    def getVerbose(self):
        """
        Return the setting of the verbosity flag.
        """
        return self.verbose
