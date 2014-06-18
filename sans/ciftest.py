#!/usr/bin/python -u
#
import sys
import os

sys.path.append( os.path.split( __file__ )[0] )

from sans import STARLexer
from sans import ErrorHandler, ContentHandler
from sans import cifparser
from sans import quote

class Test( ContentHandler, ErrorHandler ) :
    def comment( self, line, text ) :
        print "Comment:", text, "in line", line
        return False
    def startData( self, line, name ) :
        print "Start data block", name, "in line", line
        return False
    def endData( self, line, name ) :
        print "End data block", name, "in line", line
    def startSaveFrame( self, line, name ) :
        print "Start saveframe", name, "in line", line
        return False
    def endSaveFrame( self, line, name ) :
        print "End saveframe", name, "in line", line
        return False
    def startLoop( self, line ) :
        print "Start loop in line", line
        return False
    def endLoop( self, line ) :
        print "End loop in line", line
        return False
    def data( self, tag, tagline, val, valline, delim, inloop ) :
        if inloop : print "Loop",
        else : print "Free",
        print "tag/value:", tag, ":", quote( val ), "(", tagline, ":", valline, ") d", delim
        return False
    def error( self, line, msg ) :
        print "parse error in line", line, ":", msg
        return True
#
#
#
if __name__ == "__main__" :
    l = STARLexer( sys.stdin )
    t = Test()
    p = cifparser( l, t, t )
    p.parse()
#
#