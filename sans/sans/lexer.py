#!/usr/bin/python -u
#
import sys
import re

class STARLexer( object ) :
    """Regexp-based STAR lexer"""

# lexical states
    YYINITIAL = 0
    YYSEMI = 3
# tokens
    GLOBALSTART = 1
    GLOBALEND = 2
    DATASTART = 3
    DATAEND = 4
    SAVESTART = 5
    SAVEEND = 6
    LOOPSTART = 7
    STOP = 8
    TAGNAME = 9
    DVNSINGLE = 10
    DVNDOUBLE = 11
    DVNSEMICOLON = 12
    DVNFRAMECODE = 13
    DVNNON = 14
    COMMENT = 15
    FILEEND = 16
    ERROR = -1
    WARNING = -3

    _lineno = 0
    _filename = ""
    _in = None
    _buffer = ""
    _text = ""
    _yystate = YYINITIAL

    _verbose = False

    _re_comment = re.compile( "^#.*$", re.IGNORECASE )
    _re_global = re.compile( r"^global_\s*$", re.IGNORECASE )
    _re_data = re.compile( r"^data_(\S+)", re.IGNORECASE )
    _re_savestart = re.compile( r"^save_(\S+)", re.IGNORECASE )
    _re_saveend = re.compile( r"^save_", re.IGNORECASE )
    _re_loop = re.compile( "^loop_\s*$", re.IGNORECASE )
    _re_stop = re.compile( "^stop_", re.IGNORECASE )
    _re_tag = re.compile( "^_[_A-Za-z0-9]+[_.A-Za-z0-9%\-\]\[]+" )

    _re_bareword = re.compile( r"^[^\s\"';][^\s]*" )
    _re_odquote = re.compile( "^\"" )
    _re_cdquote = re.compile( "\"(?=(\\s|$))" )
    _re_osquote = re.compile( "^'" )
    _re_csquote = re.compile( r"'(?=(\s|$))" )
    _re_osemi = re.compile( "^;" )

    _re_none = re.compile( r"^\s*$" )    

    def __init__( self, inf ) :
        self._in = inf

    def _get_verbose( self ) :
        """debug flag"""
        return self._verbose
    def _set_verbose( self, flag ) :
        self._verbose = (flag == True and True or False)
    verbose = property( _get_verbose, _set_verbose )

    def getLine( self ) :
        return self._lineno

    def getText( self ) :
        return self._text

    def pushBack( self, text ) :
        if self._verbose : print "** push back |" + text + "|, buffer=|" + self._buffer + "|"
        self._buffer = text + " " + self._buffer
        if self._verbose : print "** buffer = |" + self._buffer + "|"

    def yylex( self ) :
        if self._in == None :
            print "Input file not open"
            sys.exit( 1 )
        self._text = ""

        if self._verbose : print "1) buffer: |" + self._buffer + "| len:", len( self._buffer )
        if self._buffer.isspace() : self._buffer = ""

        while len( self._buffer ) < 1 :
            line = self._in.readline()
            if self._verbose : print "line:", line
            if len( line ) == 0 : return self.FILEEND
            self._lineno = self._lineno + 1
            self._buffer = line.strip()
            if self._buffer.isspace() : self._buffer = ""

        if self._verbose : print "2) buffer: |" + self._buffer + "| len:", len( self._buffer )    

        while len( self._buffer ) > 0 :
            if self._yystate == self.YYINITIAL :
                self._buffer = self._buffer.strip()
                if len( self._buffer ) < 1 : break     # continue
# comment
                m = self._re_comment.search( self._buffer )
                if m :
                    self._text = self._buffer
                    self._buffer = ""
                    return self.COMMENT
# global
                m = self._re_global.search( self._buffer )
                if m :
                    self._buffer = self._buffer[7:]
                    self._buffer = self._buffer.strip()
                    return self.GLOBALSTART
# data
                m = self._re_data.search( self._buffer )
                if m :
                    self._text = m.group( 1 )
                    self._buffer = self._buffer[len( m.group( 0 ) ):]
                    self._buffer = self._buffer.strip()
                    return self.DATASTART
# saveframe
                m = self._re_savestart.search( self._buffer )
                if m :
                    if self._verbose : print "-- savestart", m.group( 1 )
                    self._text = m.group( 1 )
                    self._buffer = self._buffer[len( m.group( 0 ) ):]
                    self._buffer = self._buffer.strip()
                    return self.SAVESTART
                m = self._re_saveend.search( self._buffer )
                if m :
                    self._buffer = ""
                    return self.SAVEEND
# loop
                m = self._re_loop.search( self._buffer )
                if m :
                    self._buffer = self._buffer[5:]
                    self._buffer = self._buffer.strip()
                    return self.LOOPSTART
                m = self._re_stop.search( self._buffer )
                if m :
                    self._buffer = self._buffer[5:]
                    self._buffer = self._buffer.strip()
                    return self.STOP
# tag
                m = self._re_tag.search( self._buffer )
                if m :
                    self._text = m.group( 0 )
                    self._buffer = self._buffer[len( m.group( 0 ) ):]
                    self._buffer = self._buffer.strip()
                    return self.TAGNAME
# values
# dquote
                m = self._re_odquote.search( self._buffer )
                if m :
                    self._buffer = self._buffer[1:]
                    n = self._re_cdquote.search( self._buffer )
                    if not n :
                        self._text = "Unterminated double quote"
                        return self.ERROR
                    self._text = self._buffer[:n.start()]
                    self._buffer = self._buffer[n.end():]
                    self._buffer = self._buffer.strip()
                    return self.DVNDOUBLE
# squote
                m = self._re_osquote.search( self._buffer )
                if m :
                    self._buffer = self._buffer[1:]
                    n = self._re_csquote.search( self._buffer )
                    if not n :
                        self._text = "Unterminated single quote"
                        return self.ERROR
                    self._text = self._buffer[:n.start()]
                    self._buffer = self._buffer[n.end():]
                    self._buffer = self._buffer.strip()
                    return self.DVNSINGLE
# unquoted values like H' should work here
                m = self._re_bareword.search( self._buffer )
                if m :
                    self._buffer = self._buffer[len( m.group( 0 ) ):]
                    self._buffer = self._buffer.strip()
                    self._text = m.group( 0 )
                    if self._text[0] == '$' :
                        self._text = self._text[1:]
                        return self.DVNFRAMECODE
                    return self.DVNNON
# semicolon
                m = self._re_osemi.search( self._buffer )
                if m :
                    self._yystate = self.YYSEMI
                    self._buffer = self._buffer[1:]
                    if self._verbose :
                        print "entering YYSEMI, buffer:", self._buffer, "len =", len( self._buffer )
                    if not self._buffer.endswith( "\n" ) :
                        self._buffer = self._buffer + "\n" # re-add if it was stripped
                    continue
# blank
                m = self._re_none.search( self._buffer )
                if m : continue
# fail
                self._text = self._buffer 
                self._buffer = ""
                return self.ERROR

            elif self._yystate == self.YYSEMI :
                self._text = self._text + self._buffer
                if self._verbose : print "in YYSEMI, buffer:", self._buffer
# grrr. a post-condition loop would be nice. and a check for eof, too.
                while True :
                    line = self._in.readline()
                    if self._verbose : print "in YYSEMI, line:", line
                    if len( line ) == 0 : return self.FILEEND
                    self._lineno = self._lineno + 1
                    self._buffer = line              # .strip() keep original formatting
                    m = self._re_osemi.search( self._buffer )
                    if m :
                        self._buffer = self._buffer[2:]
                        if self._verbose : print "in YYSEMI, buffer:", self._buffer, "return to YYINITIAL"
                        self._yystate = self.YYINITIAL
                        return self.DVNSEMICOLON
                    else :
                        self._text = self._text + self._buffer
                        self._buffer = ""

#
#
#
def main() :
    l = STARLexer( sys.stdin )
    rc = None
    while rc != l.FILEEND :
        rc = l.yylex()
        print "RC =", rc, ":", l.getText()
#
#
#
if __name__ == "__main__" :
    main()
#
