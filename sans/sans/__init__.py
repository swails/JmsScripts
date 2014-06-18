#!/usr/bin/python -u

from __future__ import absolute_import

# suggested by one of the PEPs, probably doesn't work
#
if __name__ == "__main__" and __package__ == None :
    __package__ = "sans"

from .lexer import STARLexer
from .handlers import ErrorHandler, ContentHandler, ContentHandler2
from .SansParser import parser, parser2
from .CifParser import parser as cifparser

LONG_VALUE = 80
DEFAULT_QUOTE = STARLexer.DVNSINGLE

#
# Quote string for STAR.
#
def quote_style( value, verbose = False ) :

    import re

    global LONG_VALUE
    global DEFAULT_QUOTE

    if value == None : return STARLexer.DVNNON
    string = str( value )
    if len( string ) < 1 : return STARLexer.DVNNON

    if len( string ) > LONG_VALUE : return STARLexer.DVNSEMICOLON
    if "\n" in string : return STARLexer.DVNSEMICOLON

    dq1 = re.compile( "(^\")|(\s+\")" )
    dq2 = re.compile( "\"\s+" )
    has_dq = False
    m = dq1.search( string )
    if m : has_dq = True
    else :
        m = dq2.search( string )
        if m : has_dq = True

    sq1 = re.compile( "(^')|(\s+')" )
    sq2 = re.compile( "'\s+" )
    has_sq = False
    m = sq1.search( string )
    if m : has_sq = True
    else :
        m = sq2.search( string )
        if m : has_sq = True

    if has_sq and has_dq : return STARLexer.DVNSEMICOLON

# whitespace is only preserved inside semicolons
    string = string.strip()
    if len( string ) < 1 : return STARLexer.DVNNON

    if has_sq : return STARLexer.DVNDOUBLE
    if has_dq : return STARLexer.DVNSINGLE

    m = STARLexer._re_osemi.search( string )
    if m : return DEFAULT_QUOTE

    spc = re.compile( r"\s+" )
    m = spc.search( string )
    if m : return DEFAULT_QUOTE

    for i in ( STARLexer._re_comment, STARLexer._re_global, STARLexer._re_data, 
               STARLexer._re_saveend, STARLexer._re_loop, STARLexer._re_stop, STARLexer._re_tag ) :
        m = i.search( string )
        if m : return DEFAULT_QUOTE

    return STARLexer.DVNNON

def quote( value, verbose = False ) :

    import cStringIO

    if value == None : return "."

    qs = quote_style( value, verbose = verbose )
    if qs == STARLexer.DVNNON :
        rc = str( value ).strip()
        if len( rc ) < 1 : return "."
        else : return rc

    buf = cStringIO.StringIO()
    if qs == STARLexer.DVNSEMICOLON :
# expn: if the newline around ; was added by our pretty-printer, then it'll be adding an extra
# blank line every time the file is printed.
        if "\n" in value :
            if value.startswith( "\n" ) : buf.write( ";" )
            else : buf.write( ";\n" )
            if value.endswith( "\n" ) : buf.write( value )
            else :
                buf.write( value )
                buf.write( "\n" )
            buf.write( ";" )
        else :
            buf.write( ";\n" )
            buf.write( value )
            buf.write( "\n;" )
    elif qs == STARLexer.DVNDOUBLE :
        buf.write( '"' )
        buf.write( value )
        buf.write( '"' )
    elif qs == STARLexer.DVNSINGLE :
        buf.write( "'" )
        buf.write( value )
        buf.write( "'" )
    else :
        raise( "This can never happen" )

    rc = buf.getvalue()
    buf.close()
    return rc

__all__ = ["LONG_VALUE", "quote_style", "quote", "STARLexer", 
           "ErrorHandler", "ContentHandler", "ContentHandler2", 
           "parser", "parser2", "cifparser" ]
#
#
