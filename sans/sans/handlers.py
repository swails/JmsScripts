#!/usr/bin/python -u
#
# Interfaces for sans parsers.
#
# All callbacks can return True to stop the parser.
# File pointer remains at the start of next token when that happens.
# Callbacks must return False to keep the parser going.
#
# Content handler ones simply raise "Abstract method called"
# exception (this is good enough).
#
class ErrorHandler :
    def fatalError( self, line, msg ) :
        print "critical parse error in line", line, ":", msg
    def error( self, line, msg ) :
        print "parse error in line", line, ":", msg
        return True
    def warning( self, line, msg ) :
        print "parser warning in line", line, ":", msg
        return True

#
# This content handler returns tag/value pair in one callback
# (for both loop and free tags)
#
class ContentHandler :
    def startData( self, line, name ) :
        raise Exception( "Abstract method called" )
    def endData( self, line, name ) :
        raise Exception( "Abstract method called" )
    def startSaveframe( self, line, name ) :
        raise Exception( "Abstract method called" )
    def endSaveframe( self, line, name ) :
        raise Exception( "Abstract method called" )
    def startLoop( self, line ) :
        raise Exception( "Abstract method called" )
    def endLoop( self, line ) :
        raise Exception( "Abstract method called" )
    def comment( self, line, text ) :
        raise Exception( "Abstract method called" )
    def data( self, tag, tagline, val, valline, delim, inloop ) :
        raise Exception( "Abstract method called" )

#
# This content handler returns tag & value in separate callbacks
#
class ContentHandler2 :
    def startData( self, line, name ) :
        raise Exception( "Abstract method called" )
    def endData( self, line, name ) :
        raise Exception( "Abstract method called" )
    def startSaveframe( self, line, name ) :
        raise Exception( "Abstract method called" )
    def endSaveframe( self, line, name ) :
        raise Exception( "Abstract method called" )
    def startLoop( self, line ) :
        raise Exception( "Abstract method called" )
    def endLoop( self, line ) :
        raise Exception( "Abstract method called" )
    def comment( self, line, text ) :
        raise Exception( "Abstract method called" )
    def tag( self, line, tag ) :
        raise Exception( "Abstract method called" )
    def value( self, line, val, delim ) :
        raise Exception( "Abstract method called" )
