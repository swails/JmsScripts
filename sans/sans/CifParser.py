#!/usr/bin/python -u
#
import sys
import os
sys.path.append( os.path.split( __file__ )[0] )
from lexer import STARLexer
from handlers import ErrorHandler, ContentHandler

class parser :
    _ch = None
    _eh = None
    _lex = None
    _blockId = ""
    
    _verbose = False
    
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
            if self._verbose : print tok, self._lex.getText()
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
                self._eh.fatalError( self._lex.getLine(), "Invalid token: " + self._lex.getText() )
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
        tag = None
        tagline = -1
        val = None
        needvalue = False
        while True :
            tok = self._lex.yylex()
            if self._verbose : print "B:", tok, self._lex.getText()
            if tok == STARLexer.ERROR :
                print "+ error", self._lex.getText()
                self._eh.fatalError( self._lex.getLine(), self._lex.getText() )
                return True
            elif tok == STARLexer.WARNING :
                if self._eh.warning( self._lex.getLine(), self._lex.getText() ) :
                    return True
            elif tok == STARLexer.FILEEND :
                if needvalue :
                    if self._eh.error( self._lex.getLine(), "Premature end of file, expected value" ) :
                        return True
                self._ch.endData( self._lex.getLine(), self._blockId )
                return True
            elif tok == STARLexer.COMMENT :
                if self._ch.comment( self._lex.getLine(), self._lex.getText() ) :
                    return True
            elif tok == STARLexer.DATASTART :
                if needvalue :
                    if self._eh.error( self._lex.getLine(), "Value expected, found \"data_\"" + self._lex.getText() ) :
                        return True
# fake end of data block
                self._ch.endData( self._lex.getLine(), self._blockId )
                self._lex.pushBack( "data_" + self._lex.getText() )
                return False
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
                    if self._eh.error( self._lex.getLine(), "Value expected, found " + self._lex.getText() ) :
                        return True
                tag = self._lex.getText()
                tagline = self._lex.getLine()
                needvalue = True
            elif (tok == STARLexer.DVNSINGLE) or (tok == STARLexer.DVNDOUBLE) \
              or (tok == STARLexer.DVNNON) or (tok == STARLexer.DVNSEMICOLON) \
              or (tok == STARLexer.DVNFRAMECODE) :
                if not needvalue :
                    if self._eh.error( self._lex.getLine(), "Value not expected here: " + self._lex.getText() ) :
                        return True
                needvalue = False
                val = self._lex.getText()
                if tok == STARLexer.DVNSEMICOLON :
                    if val.startswith( "\n" ) : val = val.lstrip()

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
                if (tok == STARLexer.SAVESTART) or (tok == STARLexer.SAVEEND) :
            	    self._eh.fatalError( self._lex.getLine(), "Saveframes are not allowed in mmCIF: %s" % (self._lex.getText()) )
            	else :
            	    self._eh.fatalError( self._lex.getLine(), "Invalid token: |%s|" % (self._lex.getText()) )
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
        parsingvalues = False
        while True :
            tok = self._lex.yylex()
            if self._verbose : print "L:", tok, self._lex.getText()
            if tok == STARLexer.ERROR :
                self._eh.fatalError( self._lex.getLine(), self._lex.getText() )
                return True
            elif tok == STARLexer.WARNING :
                if self._eh.warning( self._lex.getLine(), self._lex.getText() ) :
                    return
            elif tok == STARLexer.COMMENT :
                if self._ch.comment( self._lex.getLine(), self._lex.getText() ) :
                    return True
# 
            elif (tok == STARLexer.STOP) or (tok == STARLexer.FILEEND) :
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
            	if tok == STARLexer.FILEEND :
            	    self._ch.endData( self._lex.getLine(), self._blockId )
            	    return True
                return False
# real end of loop variants: data_, loop_, or tag
# -- return to data block level
            elif tok == STARLexer.DATASTART :
                if parsingtags :
                    if self._eh.error( self._lex.getLine(), "Tag or value expected, found \"data_\"" + self._lex.getText() ) :
                        return True
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
# fake end of loop
                if self._ch.endLoop( self._lex.getLine() ) :
                    return True
# don't fake end of data block -- let parseDataBlock() do it
#                self._ch.endData( self._lex.getLine(), self._blockId )
                self._lex.pushBack( "data_" + self._lex.getText() )
                return False

            elif tok == STARLexer.LOOPSTART :
                if parsingtags :
                    if self._eh.error( self._lex.getLine(), "Tag or value expected, found \"loop_\"" ) :
                        return True
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
# fake end of loop
                if self._ch.endLoop( self._lex.getLine() ) :
                    return True
                self._lex.pushBack( "loop_" )
                return False

            elif tok == STARLexer.TAGNAME :
                if parsingvalues :
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
# fake end of loop
                    if self._ch.endLoop( self._lex.getLine() ) :
                        return True
                    self._lex.pushBack( self._lex.getText() )
                    return False
# else add to the list
                tags.append( self._lex.getText() )
                taglines.append( self._lex.getLine() )

            elif (tok == STARLexer.DVNSINGLE) or (tok == STARLexer.DVNDOUBLE) \
              or (tok == STARLexer.DVNNON) or (tok == STARLexer.DVNSEMICOLON) \
              or (tok == STARLexer.DVNFRAMECODE) :
                if len( tags ) < 1 :
                    if self._eh.error( self._lex.getLine(), "Loop with no tags: expected tag, found" + self._lex.getText() ) :
                        return True
                if parsingtags : parsingtags = False
                parsingvalues = True
                val = self._lex.getText()
                if tok == STARLexer.DVNSEMICOLON :
                    if val.startswith( "\n" ) : val = val.lstrip()
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
                self._eh.fatalError( self._lex.getLine(), "Invalid token: " + self._lex.getText() )
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
