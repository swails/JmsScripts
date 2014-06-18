#!/usr/bin/python -u
#
import sys
import os
sys.path.append( os.path.split( __file__ )[0] )
from lexer import STARLexer
from handlers import ErrorHandler, ContentHandler

#
# returns tag/value (loop or free) pair in one callback
#    
class parser :
    _ch = None
    _eh = None
    _lex = None
    _blockId = ""

    def __init__( self, lex, ch, eh ) :
        self._lex = lex
        self._ch = ch
        self._eh = eh

    def parse( self ) :
        if self._lex == None : 
            print "Lexer not initialized"
            sys.exit( 1 )
        if self._ch == None : 
            print "Content handler not initialized"
            sys.exit( 1 )
        if self._eh == None : 
            print "Error handler not initialized"
            sys.exit( 1 )
        while True :
            tok = self._lex.yylex()
            if tok == STARLexer.ERROR :
                self._eh.fatalError( self._lex.getLine(), self._lex.getText() )
                return
            elif tok == STARLexer.WARNING :
                if self._eh.warning( self._lex.getLine(), self._lex.getText() ) :
                    return
            elif tok == STARLexer.FILEEND :
                self._ch.endData( self._lex.getLine(), self._blockId )
                return
            elif tok == STARLexer.COMMENT :
                if self._ch.comment( self._lex.getLine(), self._lex.getText() ) :
                    return
            elif tok == STARLexer.DATASTART :
                self._blockId = self._lex.getText()
                if self._ch.startData( self._lex.getLine(), self._blockId ) :
                    return
                if self.parseDataBlock() : 
                    return
            else :
                self._eh.fatalError( self._lex.getLine(), "Invalid token: %s" % (self._lex.getText()) )
                return

#
#
#
    def parseDataBlock( self ) :
        if self._lex == None : 
            print "Lexer not initialized"
            sys.exit( 1 )
        if self._ch == None : 
            print "Content handler not initialized"
            sys.exit( 1 )
        if self._eh == None : 
            print "Error handler not initialized"
            sys.exit( 1 )
        while True :
            tok = self._lex.yylex()
            if tok == STARLexer.ERROR :
                print "+ error", self._lex.getText()
                self._eh.fatalError( self._lex.getLine(), self._lex.getText() )
                return True
            elif tok == STARLexer.WARNING :
                if self._eh.warning( self._lex.getLine(), self._lex.getText() ) :
                    return True
            elif tok == STARLexer.FILEEND :
                self._ch.endData( self._lex.getLine(), self._blockId )
                return True
            elif tok == STARLexer.COMMENT :
                if self._ch.comment( self._lex.getLine(), self._lex.getText() ) :
                    return True
            elif tok == STARLexer.SAVESTART :
                if self._ch.startSaveFrame( self._lex.getLine(), self._lex.getText() ) :
                    return True
                if self.parseSaveFrame() : 
                    return True
            else :
                self._eh.fatalError( self._lex.getLine(), "Invalid token in data block: %s" % (self._lex.getText()) )
                return True

#
#
#
    def parseSaveFrame( self ) :
        if self._lex == None : 
            print "Lexer not initialized"
            sys.exit( 1 )
        if self._ch == None : 
            print "Content handler not initialized"
            sys.exit( 1 )
        if self._eh == None : 
            print "Error handler not initialized"
            sys.exit( 1 )
        tag = None
        tagline = -1
        val = None
        needvalue = False
        while True :
            tok = self._lex.yylex()
            if tok == STARLexer.ERROR :
                self._eh.fatalError( self._lex.getLine(), self._lex.getText() )
                return True
            elif tok == STARLexer.WARNING :
                if self._eh.warning( self._lex.getLine(), self._lex.getText() ) :
                    return
            elif tok == STARLexer.FILEEND :
                if self._eh.error( self._lex.getLine(), "Premature end of file: no closing \"save_\"" ) :
                    return True
                self._ch.endData( self._lex.getLine(), self._blockId )
                return True
            elif tok == STARLexer.STOP :
                self._eh.fatalError( self._lex.getLine(), "Found \"stop_\", expected \"save_\"" )
                return True
            elif tok == STARLexer.COMMENT :
                if self._ch.comment( self._lex.getLine(), self._lex.getText() ) :
                    return True
            elif tok == STARLexer.SAVEEND :
                if needvalue :
                    if self._eh.error( self._lex.getLine(), "Value expected, found \"save_\"" ) :
                        return True
                return self._ch.endSaveFrame( self._lex.getLine(), self._lex.getText() )
            elif tok == STARLexer.LOOPSTART :
                if needvalue :
                    if self._eh.error( self._lex.getLine(), "Value expected, found \"loop_\"" ) :
                        return True
                if self._ch.startLoop( self._lex.getLine() ) :
                    return True
                if self.parseLoop() :
                    return True
            elif tok == STARLexer.TAGNAME :
                if needvalue :
                    if self._eh.error( self._lex.getLine(), "Value expected, found %s" % (self._lex.getText()) ) :
                        return True
                tag = self._lex.getText()
                tagline = self._lex.getLine()
                needvalue = True
            elif tok in (STARLexer.DVNSINGLE, STARLexer.DVNDOUBLE, STARLexer.DVNNON, STARLexer.DVNSEMICOLON, STARLexer.DVNFRAMECODE) :
                if not needvalue :
                    if self._eh.error( self._lex.getLine(), "Value not expected here: %s" % (self._lex.getText()) ) :
                        return True
                needvalue = False
                val = self._lex.getText()
                if tok in ( STARLexer.DVNSINGLE, STARLexer.DVNDOUBLE ) :
                    if ((val[0] in "'") and (val[-1] == "'") or (val[0] in '"') and (val[-1] == '"')) :
                        if self._eh.warning( self._lex.getLine(), "value already in quotes: %s" % (val)  ) :
                            return True
                if tok == STARLexer.DVNSEMICOLON :
                    if val.startswith( "\n" ) : val = val.lstrip( "\n" )
                    n = STARLexer._re_data.search( val )
                    if not n : n = STARLexer._re_saveend.search( val )
                    if not n : n = STARLexer._re_loop.search( val )
                    if not n : n = STARLexer._re_stop.search( val )
                    if not n : n = STARLexer._re_tag.search( val )
                    if n :
                        if self._eh.warning( self._lex.getLine(), "Keyword in value: %s" % (n.group()) ) :
                            return True
                if len( val ) < 1 : 
                    val = None
                    if self._eh.error( self._lex.getLine(), "NULL value!" ) :
                        return True

                if self._ch.data( tag, tagline, val, self._lex.getLine(), tok, False ) :
                    return True
                tag = None
                tagline = -1
                val = None
            else :
                self._eh.fatalError( self._lex.getLine(), "Invalid token in saveframe: %s" % (self._lex.getText()) )
                return True

#
#
#
    def parseLoop( self ) :
        if self._lex == None : 
            print "Lexer not initialized"
            sys.exit( 1 )
        if self._ch == None : 
            print "Content handler not initialized"
            sys.exit( 1 )
        if self._eh == None : 
            print "Error handler not initialized"
            sys.exit( 1 )
        tags = []
        taglines = []
        numvals = 0
        loopcol = 0
        lastline = -1
        wrongline = -1
        wrongcol = -1
        val = None
        tag = None
        tagline = -1
        parsingtags = True
        while True :
            tok = self._lex.yylex()
            if tok == STARLexer.ERROR :
                self._eh.fatalError( self._lex.getLine(), self._lex.getText() )
                return True
            elif tok == STARLexer.WARNING :
                if self._eh.warning( self._lex.getLine(), self._lex.getText() ) :
                    return
            elif tok == STARLexer.FILEEND :
                if self._eh.error( self._lex.getLine(), "Premature end of file: no closing \"stop_\"" ) :
                    return True
                self._ch.endData( self._lex.getLine(), self._blockId )
                return True
            elif tok == STARLexer.SAVEEND :
                self._eh.fatalError( self._lex.getLine(), "Found \"save_\", expected \"stop_\"" )
                return True
            elif tok == STARLexer.LOOPSTART :
                self._eh.fatalError( self._lex.getLine(), "\"loop_\" not expected here" )
                return True
            elif tok == STARLexer.COMMENT :
                if self._ch.comment( self._lex.getLine(), self._lex.getText() ) :
                    return True
# normal exit point
            elif tok == STARLexer.STOP :
# loop checks
                if len( tags ) < 1 :
                    if self._eh.error( self._lex.getLine(), "Loop with no tags" ) :
                        return True
                if numvals < 1 :
                    if self._eh.error( self._lex.getLine(), "Loop with no values" ) :
                        return True
# loop count
                if (numvals % len( tags )) != 0 :
                    if wrongline < 0 : wrongline = self._lex.getLine()
                    if self._eh.warning( wrongline, "Loop count error" ) :
                        return True
                if self._ch.endLoop( self._lex.getLine() ) :
                    return True
                return False
            elif tok == STARLexer.TAGNAME :
                if not parsingtags :
                    if self._eh.error( self._lex.getLine(), "Value expected, found %s" % (self._lex.getText()) ) :
                        return True
                tags.append( self._lex.getText() )
                taglines.append( self._lex.getLine() )
            elif tok in (STARLexer.DVNSINGLE, STARLexer.DVNDOUBLE, STARLexer.DVNNON, STARLexer.DVNSEMICOLON, STARLexer.DVNFRAMECODE) :
                if len( tags ) < 1 :
                    if self._eh.error( self._lex.getLine(), "Loop with no tags: expected tag, found %s" % (self._lex.getText()) ) :
                        return True
                if parsingtags : parsingtags = False
                val = self._lex.getText()
                if tok in ( STARLexer.DVNSINGLE, STARLexer.DVNDOUBLE ) :
                    if ((val[0] in "'") and (val[-1] == "'") or (val[0] in '"') and (val[-1] == '"')) :
                        if self._eh.warning( self._lex.getLine(), "value already in quotes: %s" % (val)  ) :
                            return True
                if tok == STARLexer.DVNSEMICOLON :
                    if val.startswith( "\n" ) : val = val.lstrip( "\n" )
                    n = STARLexer._re_data.search( val )
                    if not n : n = STARLexer._re_saveend.search( val )
                    if not n : n = STARLexer._re_loop.search( val )
                    if not n : n = STARLexer._re_stop.search( val )
                    if not n : n = STARLexer._re_tag.search( val )
                    if n :
                        if self._eh.warning( self._lex.getLine(), "Keyword in value: %s" % (n.group())  ) :
                            return True
                numvals += 1
                tag = tags[loopcol]
                tagline = taglines[loopcol]
                loopcol += 1
                if loopcol == len( tags ) :
                    if lastline != self._lex.getLine() :
                        if wrongline < 0 : wrongline = self._lex.getLine()
                        lastline = self._lex.getLine()
                    loopcol = 0
                if len( val ) < 1 : 
                    val = None
                    if self._eh.error( self._lex.getLine(), "NULL value!" ) :
                        return True
                if self._ch.data( tag, tagline, val, self._lex.getLine(), tok, True ) :
                    return True
                tag = None
                tagline = -1
                val = None
            else :
                self._eh.fatalError( self._lex.getLine(), "Invalid token in loop: %s" % (self._lex.getText()) )
                return True

#
# returns tag and value in separate callbacks
# 
class parser2 :
    _ch = None
    _eh = None
    _lex = None
    _blockId = ""

    def __init__( self, lex, ch, eh ) :
        self._lex = lex
        self._ch = ch
        self._eh = eh

    def parse( self ) :
        if self._lex == None : 
            print "Lexer not initialized"
            sys.exit( 1 )
        if self._ch == None : 
            print "Content handler not initialized"
            sys.exit( 1 )
        if self._eh == None : 
            print "Error handler not initialized"
            sys.exit( 1 )
        while True :
            tok = self._lex.yylex()
            if tok == STARLexer.ERROR :
                self._eh.fatalError( self._lex.getLine(), self._lex.getText() )
                return
            elif tok == STARLexer.WARNING :
                if self._eh.warning( self._lex.getLine(), self._lex.getText() ) :
                    return
            elif tok == STARLexer.FILEEND :
                self._ch.endData( self._lex.getLine(), self._blockId )
                return
            elif tok == STARLexer.COMMENT :
                if self._ch.comment( self._lex.getLine(), self._lex.getText() ) :
                    return
            elif tok == STARLexer.DATASTART :
                self._blockId = self._lex.getText()
                if self._ch.startData( self._lex.getLine(), self._blockId ) :
                    return
                if self.parseDataBlock() : 
                    return
            else :
                self._eh.fatalError( self._lex.getLine(), "Invalid token: %s" % (self._lex.getText()) )
                return

#
#
#
    def parseDataBlock( self ) :
        if self._lex == None : 
            print "Lexer not initialized"
            sys.exit( 1 )
        if self._ch == None : 
            print "Content handler not initialized"
            sys.exit( 1 )
        if self._eh == None : 
            print "Error handler not initialized"
            sys.exit( 1 )
        while True :
            tok = self._lex.yylex()
            if tok == STARLexer.ERROR :
                print "+ error", self._lex.getText()
                self._eh.fatalError( self._lex.getLine(), self._lex.getText() )
                return True
            elif tok == STARLexer.WARNING :
                if self._eh.warning( self._lex.getLine(), self._lex.getText() ) :
                    return True
            elif tok == STARLexer.FILEEND :
                self._ch.endData( self._lex.getLine(), self._blockId )
                return True
            elif tok == STARLexer.COMMENT :
                if self._ch.comment( self._lex.getLine(), self._lex.getText() ) :
                    return True
            elif tok == STARLexer.SAVESTART :
                if self._ch.startSaveFrame( self._lex.getLine(), self._lex.getText() ) :
                    return True
                if self.parseSaveFrame() : 
                    return True
            else :
                self._eh.fatalError( self._lex.getLine(), "Invalid token in data block: %s" % (self._lex.getText()) )
                return True

#
#
#
    def parseSaveFrame( self ) :
        if self._lex == None : 
            print "Lexer not initialized"
            sys.exit( 1 )
        if self._ch == None : 
            print "Content handler not initialized"
            sys.exit( 1 )
        if self._eh == None : 
            print "Error handler not initialized"
            sys.exit( 1 )
        tag = None
        tagline = -1
        val = None
        needvalue = False
        while True :
            tok = self._lex.yylex()
            if tok == STARLexer.ERROR :
                self._eh.fatalError( self._lex.getLine(), self._lex.getText() )
                return True
            elif tok == STARLexer.WARNING :
                if self._eh.warning( self._lex.getLine(), self._lex.getText() ) :
                    return
            elif tok == STARLexer.FILEEND :
                if self._eh.error( self._lex.getLine(), "Premature end of file: no closing \"save_\"" ) :
                    return True
                self._ch.endData( self._lex.getLine(), self._blockId )
                return True
            elif tok == STARLexer.STOP :
                self._eh.fatalError( self._lex.getLine(), "Found \"stop_\", expected \"save_\"" )
                return True
            elif tok == STARLexer.COMMENT :
                if self._ch.comment( self._lex.getLine(), self._lex.getText() ) :
                    return True
            elif tok == STARLexer.SAVEEND :
                if needvalue :
                    if self._eh.error( self._lex.getLine(), "Value expected, found \"save_\"" ) :
                        return True
                return self._ch.endSaveFrame( self._lex.getLine(), self._lex.getText() )
            elif tok == STARLexer.LOOPSTART :
                if needvalue :
                    if self._eh.error( self._lex.getLine(), "Value expected, found \"loop_\"" ) :
                        return True
                if self._ch.startLoop( self._lex.getLine() ) :
                    return True
                if self.parseLoop() :
                    return True
            elif tok == STARLexer.TAGNAME :
                if needvalue :
                    if self._eh.error( self._lex.getLine(), "Value expected, found \"loop_\"" ) :
                        return True
                if self._ch.tag( self._lex.getLine(), self._lex.getText() ) :
                    return True
                needvalue = True
            elif tok in (STARLexer.DVNSINGLE, STARLexer.DVNDOUBLE, STARLexer.DVNNON, STARLexer.DVNSEMICOLON, STARLexer.DVNFRAMECODE) :
                if not needvalue :
                    if self._eh.error( self._lex.getLine(), "Value not expected here: %s" % (self._lex.getText()) ) :
                        return True
                needvalue = False
                val = self._lex.getText()
                if tok in ( STARLexer.DVNSINGLE, STARLexer.DVNDOUBLE ) :
                    if ((val[0] in "'") and (val[-1] == "'") or (val[0] in '"') and (val[-1] == '"')) :
                        if self._eh.warning( self._lex.getLine(), "value already in quotes: %s" % (val)  ) :
                            return True
                if tok == STARLexer.DVNSEMICOLON :
                    if val.startswith( "\n" ) : val = val.lstrip( "\n" )
                    n = STARLexer._re_data.search( val )
                    if not n : n = STARLexer._re_saveend.search( val )
                    if not n : n = STARLexer._re_loop.search( val )
                    if not n : n = STARLexer._re_stop.search( val )
                    if not n : n = STARLexer._re_tag.search( val )
                    if n :
                        if self._eh.warning( self._lex.getLine(), "Keyword in value: %s" % (n.group())  ) :
                            return True
                if len( val ) < 1 : 
                    val = None
                    if self._eh.error( self._lex.getLine(), "NULL value!" ) :
                        return True
                if self._ch.value( self._lex.getLine(), val, tok ) :
                    return True
                val = None
            else :
                self._eh.fatalError( self._lex.getLine(), "Invalid token in saveframe: %s" % (self._lex.getText()) )
                return True

#
#
#
    def parseLoop( self ) :
        if self._lex == None : 
            print "Lexer not initialized"
            sys.exit( 1 )
        if self._ch == None : 
            print "Content handler not initialized"
            sys.exit( 1 )
        if self._eh == None : 
            print "Error handler not initialized"
            sys.exit( 1 )
        numtags = 0
        numvals = 0
        loopcol = 0
        lastline = -1
        wrongline = -1
        wrongcol = -1
        parsingtags = True
        val = None
        while True :
            tok = self._lex.yylex()
            if tok == STARLexer.ERROR :
                self._eh.fatalError( self._lex.getLine(), self._lex.getText() )
                return True
            elif tok == STARLexer.WARNING :
                if self._eh.warning( self._lex.getLine(), self._lex.getText() ) :
                    return
            elif tok == STARLexer.FILEEND :
                if self._eh.error( self._lex.getLine(), "Premature end of file: no closing \"stop_\"" ) :
                    return True
                self._ch.endData( self._lex.getLine(), self._blockId )
                return True
            elif tok == STARLexer.SAVEEND :
                self._eh.fatalError( self._lex.getLine(), "Found \"save_\", expected \"stop_\"" )
                return True
            elif tok == STARLexer.LOOPSTART :
                self._eh.fatalError( self._lex.getLine(), "\"loop_\" not expected here" )
                return True
            elif tok == STARLexer.COMMENT :
                if self._ch.comment( self._lex.getLine(), self._lex.getText() ) :
                    return True
# normal exit point
            elif tok == STARLexer.STOP :
# loop checks
                if numtags < 1 :
                    if self._eh.error( self._lex.getLine(), "Loop with no tags" ) :
                        return True
                if numvals < 1 :
                    if self._eh.error( self._lex.getLine(), "Loop with no values" ) :
                        return True
# loop count
                if (numvals % numtags) != 0 :
                    if wrongline < 0 : wrongline = self._lex.getLine()
                    if self._eh.warning( wrongline, "Loop count error" ) :
                        return True
                if self._ch.endLoop( self._lex.getLine() ) :
                    return True
                return False
            elif tok == STARLexer.TAGNAME :
                if not parsingtags :
                    if self._eh.error( self._lex.getLine(), "Value expected, found %s" % (self._lex.getText()) ) :
                        return True
                numtags += 1
                if self._ch.tag( self._lex.getLine(), self._lex.getText() ) :
                    return True

            elif tok in (STARLexer.DVNSINGLE, STARLexer.DVNDOUBLE, STARLexer.DVNNON, STARLexer.DVNSEMICOLON, STARLexer.DVNFRAMECODE) :
                if numtags < 1 :
                    if self._eh.error( self._lex.getLine(), "Loop with no tags: expected tag, found %s" % (self._lex.getText()) ) :
                        return True
                if parsingtags : parsingtags = False
                val = self._lex.getText()
                if tok in ( STARLexer.DVNSINGLE, STARLexer.DVNDOUBLE ) :
                    if ((val[0] in "'") and (val[-1] == "'") or (val[0] in '"') and (val[-1] == '"')) :
                        if self._eh.warning( self._lex.getLine(), "value already in quotes: %s" % (val)  ) :
                            return True
                if tok == STARLexer.DVNSEMICOLON :
                    if val.startswith( "\n" ) : val = val.lstrip( "\n" )
                    n = STARLexer._re_data.search( val )
                    if not n : n = STARLexer._re_saveend.search( val )
                    if not n : n = STARLexer._re_loop.search( val )
                    if not n : n = STARLexer._re_stop.search( val )
                    if not n : n = STARLexer._re_tag.search( val )
                    if n :
                        if self._eh.warning( self._lex.getLine(), "Keyword in value: %s" % (n.group())  ) :
                            return True
                numvals += 1
                loopcol += 1
                if loopcol == numtags :
                    if lastline != self._lex.getLine() :
                        if wrongline < 0 : wrongline = self._lex.getLine()
                        lastline = self._lex.getLine()
                    loopcol = 0
                if len( val ) < 1 : 
                    val = None
                    if self._eh.error( self._lex.getLine(), "NULL value!" ) :
                        return True
                if self._ch.value( self._lex.getLine(), val, tok ) :
                    return True
                val = None
            else :
                self._eh.fatalError( self._lex.getLine(), "Invalid token in loop: %s" % (self._lex.getText()) )
                return True

#
#
#
def main() :
    l = STARLexer( sys.stdin )
    e = ErrorHandler()
    c = ContentHandler()
    p = parser( l, c, e )
    p.parse()
#
#
#
if __name__ == "__main__" :
    main()
#
#