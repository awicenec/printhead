# printhead

FITS Header tool in pure python

This small tool has a long history and come a long way from early versions of Python to the latest Python 3. Its original intend was to
allow users to display FITS headers on screen using just a single shell command. Later a lot of additional functionality was added.

## Installation

Since it is pure python without any external dependencies installation is very simple.
From the main directory where also this README file is located run:
```
pip install .
```

## Short User Guide

The short help, which is displayed when calling it without any parameters like ./printhead.py, gives a brief overview of the functionality.
```
Script dumps headers of FITS files to stdout or creates header
files. It supports compressed files (.gz and .Z)

Synopsis: printhead.py [-s <KEYWORD> -H <number> -M <extnum> -S -x <type> -e -h]
 fname1 [fname2]...

If only a file name is given, the primary header of that file is printed.

--extract|-e:   All the headers of the files found are then
                extracted to directories with the same name as the last
                directory found in the path-names of the files. The
                header files will have the same base name as the file, but
                with the extension '.hdr'.
--skey|-s:      if the given KEYWORD is found
                in the header only the matching lines will be printed.

--header|-H:    <number> specifies the number of the header to be printed.
                If 99 is given, all headers are printed. If <number> is
                negative only the structure of the file is printed.

--xml|-x:       <type> is either 'vo' or 'xf'. All the headers of the files found
                are then extracted to directories with the same name as the last
                directory found in the path-names of the files. The
                header files will have the same base name as the file, but
                with the extension '.xml'. The files use XFits as a format
                if 'xf' is specified and VOTable format if 'vo' is specified
                for <type>.
--Struct|-S     Show the structure of the FITS file

--tsv|-t        Print keywords as a tab-separated-value list. The list contains
                the items: keyword name, keyword value, comment, keyword type, index
                The index item is the running number of the keyword within the header.

--check|-c      If this flag is set the CRC32 checksum of the data part of the
                extensions is calculated.
--mode-m        [1]|0: If set to 0 the program does not try to skip the data part
                between headers. This is useful for interpreting header files.
--parse|-p      Switch the full parsing of the header on
                extensions is calculated.
--help|-h:      print this help and exit.

Version: 5.0
```
