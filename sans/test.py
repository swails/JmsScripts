#!/usr/bin/python -u
#
import sys
import os
sys.path.append( os.path.split( __file__ )[0] )

import sans

class Test( sans.ContentHandler, sans.ErrorHandler ) :
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
        print "tag/value:", tag, ":", sans.quote( val ), "(", tagline, ":", valline, ") d", delim
        return False
    def fatalError( self, line, msg ) :
        print "parse fatal error in line", line, ":", msg
        return True
    def error( self, line, msg ) :
        print "parse error in line", line, ":", msg
        return True
    def warning( self, line, msg ) :
        print "parse warning in line", line, ":", msg
        return False

#
#
#
if __name__ == "__main__" :
    l = sans.STARLexer( sys.stdin )
    t = Test()
    p = sans.parser( l, t, t )
    p.parse()
#
#
